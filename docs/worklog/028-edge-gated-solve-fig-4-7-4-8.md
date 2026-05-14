# Worklog 028 — Edge-gated pole-field solve + FW 2011 Fig 4.7 / 4.8

**Date**: 2026-05-14 (continues worklog 027)
**Author**: Claude Opus
**Beads**: `padetaylor-dmb` (edge-gated solve), `padetaylor-7na`
(Fig 4.7/4.8) — both closed.
**Scope**: Reproduce FW 2011 Fig 4.7 (six PI pole fields) and Fig 4.8
(PI pole field with a sharp pattern transition). Started as a figure
task; uncovered a genuine solver limitation that needed a new module
to fix.

> **Take-home**: `path_network_solve` is an IVP solver and FW 2011
> md:401 is blunt — "smooth regions are unstable regions for any IVP
> solver". Integrate the PI *tritronquée* across a large grid and
> accumulated error perturbs it off its (four-of-five-sectors-empty)
> manifold, blooming spurious poles where there should be none. The
> fix is FW's own: gate the IVP so it never targets a smooth cell.
> `src/EdgeGatedSolve.jl` (new) does this by region growing — expand
> the target set only into edge-detector-confirmed cells *connected*
> to the field already found. Both figures ship; test suite **1509 →
> 1521 GREEN** (+12 from `edge_gated_solve_test.jl`).

## How the figure task turned into a solver task

The plan was simple: `path_network_solve` over a grid → `extract_poles`
(worklog 027) → scatter plot. A probe of panel (f), the tritronquée,
killed that plan: over `[-50,50]²` the extracted pole field was nearly
angle-uniform — the four pole-free sectors had filled with spurious
poles. Over `[-12,12]²` about a quarter of the poles still leaked into
them.

This is not a bug to patch — it is the documented nature of the
method. FW 2011 says so directly:

  - **md:401** — "smooth regions are unstable regions for any IVP
    solver, and the associated loss of accuracy ought not be carried
    back into the pole field. One way to prevent this is to force the
    path selection algorithm to complete pole fields before stepping
    into smooth regions … with the pole field edge detection procedure
    described in Section 3.2.2, this process can readily be fully
    automated."
  - **md:147** — Fig 4.7 panels **a–e** are pure-IVP "5×5 composites
    of 25 independent runs"; **(f)** is pointedly *not* in that list.
  - **md:216-224** — the tritronquée family is handled with a BVP
    solver for the smooth band + IVP only "to run out the pole fields".
  - **md:200, 350** — FW's pole-field integrator runs in plain MATLAB
    double precision; extended precision (Maple, 32 digits) is used
    *only* to get the tritronquée's IC values.

So the cure is *architectural*, not higher precision: confine the IVP
to the pole field. A probe confirmed it — restricting the path-network
targets to the (known) tritronquée wedge gave a pole field with the
four other sectors *exactly* empty.

## `src/EdgeGatedSolve.jl` (new, ~110 code LOC + literate docstring)

`edge_gated_pole_field_solve(prob, xs, ys; …) -> EdgeGatedSolution` —
region growing:

  1. **Seed** `field` ← grid cells within `seed_radius` of the IC.
  2. **Grow** to a fixpoint: dilate `field` by `grow_rings`,
     `path_network_solve` the dilated target set, `pole_field_mask`
     (FW eq. 3.3), **clean the mask**, admit newly-reached cells.
  3. **Final solve** over `field` + one thin ring; its
     `PathNetworkSolution` feeds `extract_poles`.

### The non-obvious part — the edge mask alone is not a usable gate

FW's `level = 0.001` threshold is calibrated for a *fine* grid (FW
used 161×161 over `[-10,10]²`, spacing ≈ 0.125). The smooth-region
discretisation residual of the 5-point stencil scales as
`Δu ∼ h²·u⁗/12`, so on a coarse lattice it is orders of magnitude
larger and the per-cell classifier is *noisy*: probes at spacing 2
showed the wedge and smooth `log₁₀|Δu|` distributions essentially
overlapping. Two cheap morphological steps turn the noisy classifier
into a reliable region gate:

  - **Opening** (erode then dilate) removes false-positive specks and
    one-cell-wide bridges; the genuine pole field is a thick connected
    blob and survives.
  - **Flood-fill** from the current `field` keeps only mask cells
    *connected* to it — a stranded false positive, or a thin bridge
    the opening missed, is never reached, so the field cannot leak
    across a smooth gap.

This pair is the automatic stand-in for FW's "we implemented this
manually after having inspected … where the pole field edges
appeared" (md:401). It still needs a lattice spacing `≲ 1` — at
spacing ≈ 2 even opening + flood-fill cannot recover a clean gate, and
the module docstring says so ("Grid resolution matters") rather than
silently accepting a doomed grid.

### TDD record (`test/edge_gated_solve_test.jl`, 12 assertions)

Oracle is the PI tritronquée. EG.1.1 (the solve grew a region),
EG.1.2 (sector confinement — the four pole-free sectors stay
near-empty: `nfree ≤ 5`, `nfree < 0.1·npop`), EG.1.3 (**the
load-bearing contrast** — the plain ungated solve over the same grid
*does* leak: `nfree_gated < nfree_plain ÷ 4`), EG.2.1 (input
validation). Mutation-proven — M1 (drop the connected-mask gate), M2
(disable region growing), M3 (thicken the final frontier ring); each
confirmed RED, reverted. Procedure in the test file footer.

## The figures

`figures/fw2011_fig_4_7.jl` (6 panels, edge-gated, 101×101 over
`[-50,50]²`) and `figures/fw2011_fig_4_8.jl` (edge-gated, 121×61 over
`[-90,30]×[-30,30]`). A `pole_scatter_axis` helper was added to
`figures/figutil.jl`.

All six Fig 4.7 panels use `edge_gated_pole_field_solve` uniformly:

  - **(a, c, d)** — generic solutions, poles in every sector. The
    edge gate classifies essentially the whole window as pole-field
    (≈ 96% of cells) and the region grows to fill it — i.e. the gated
    solve reduces to the plain pure-IVP fill FW used for a–e, reached
    safely. Visual match: (a)'s sparse diagonal seam, (d)'s uniform
    dense field.
  - **(b)** — near-tronquée: the gated solve resolves the curved
    smooth band sweeping from the upper-left, matching FW's (b).
  - **(e)** — near-tritronquée: the X-shaped structure with smooth
    wedges near the origin, matching FW's (e).
  - **(f)** — tritronquée: the gate confines the solve to **1847 of
    10201 cells** (≈ 18 %, one of five sectors). The extracted pole
    field is a clean wedge around the positive real axis with the
    other four sectors *empty* — the angular histogram is
    `[0,0,0,0,245,929,933,232,0,0,0,0]`. This is the panel a plain
    solve gets qualitatively wrong, and the reason the module exists.

`figures/fw2011_fig_4_8.jl` reproduces the `u(0)=-5, u'(0)=0` pole
field; the on-axis poles bracket FW's stated transition at `Re ≈ -60`
(nearest on-axis poles at −63.9 and −59.7, with the pattern texture
visibly changing across that band).

Both are *visual* reproductions, the same acceptance bar Fig 4.3/4.4
shipped under (worklog 024); the quantitative per-panel pole-count pin
is left to bead `padetaylor-p3l`'s test-side scope.

## Frictions / deferred

  - **BVP is not composed here.** The user's steer was "compose
    EdgeDetector + BVP". For the tritronquée the smooth sectors are
    *unbounded* — not the horizontally-BVP-bridgeable geometry
    `lattice_dispatch_solve` handles (FW md:190), and for a
    pole-*location* plot the smooth sectors carry no poles, so leaving
    them empty is correct. BVP-fill of *bounded* smooth bands stays
    with `lattice_dispatch_solve`; the tritronquée's imaginary-axis
    BVP is the separate Fig 4.1 pattern, bead `padetaylor-gky`. This
    is recorded honestly in `EdgeGatedSolve`'s docstring rather than
    forcing a tool where it does not fit (Rule 9).
  - **Region growing re-solves from scratch each pass.** `path_network_
    solve` has no warm start, so each growth pass rebuilds the tree.
    For the figure scripts (run once) the ~90–110 s/panel cost is
    acceptable; a warm-startable path-network is a possible future
    optimisation, not filed.

## Hard-won lesson

**A "reproduce the figure" task can be load-bearing on a capability
that does not exist yet.** The honest read of FW md:401 — *smooth
regions are unstable for any IVP solver* — is that the plain
path-network simply cannot produce Fig 4.7(f), and no amount of
tolerance-fiddling changes that. The fix had to be the architectural
one FW themselves describe. Reading the paper's *limitations* as
carefully as its method (Law 1) is what turned a stuck figure into a
shipped module.
