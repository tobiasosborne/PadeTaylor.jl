# Worklog 063 — zoom-figure iteration arc + Stage-2 fill diagnosis

**Date**: 2026-05-27
**Author**: Claude Opus 4.7
**Epic**: `padetaylor-0ln` (v0.2) · **Follows**: worklog 062
**Scope**: Iterating the two complex-function zoom figures introduced in
worklog 062 (`figures/kkg_pi2_abs_heatmap_zoom.{jl,png}` and
`figures/kkg_pi2_abs_phase_surface_zoom.{jl,png}`) — first from `[-2, 2]²`
@ N=801 (user-approved) to `[-3, 3]²` @ N=1001 (user-requested), then
through four iterations on Stage-2 fill quality. Final ship: 21×21=441
dense target lattice, walk step `h = 0.05`, `extrapolate = false` +
BFS nearest-non-NaN fill in the rendering scripts.

> **Take-home**: zoom over `[-3, 3]²` necessarily includes inner-wedge
> pole territory `|x| ∈ [2, 3], |arg x| < 36°`. The zoom kernel uses a
> single Padé walk + single-nearest-node Stage-2 fill (no VC-4 cross-
> validation, deliberate simplification from the parent surface kernel)
> — that fill produces some artefacts near pole regions where raw Padé
> Q-roots stand in for real poles. The four iterations documented here
> minimised but did not eliminate those artefacts. The principled fix
> is K-nearest cross-validated Stage-2 fill in `src/VectorPathNetwork.jl`,
> deferred as a follow-up.

## The iteration arc

### Iteration 1 — `[-3, 3]²` @ N=1001, 8 targets, `extrapolate = true`

Started by changing `ZOOM_XY_LIM = 2.0 → 3.0` and `ZOOM_GRID_N = 801 →
1001` in `figures/_kkg_pi2_zoom_helpers.jl`. Kept the 8-perimeter target
list (corners + edge midpoints) the [-2, 2]² version used.

**Result**: kernel runtime 60.4 s (cache MISS, fresh compute), 221
visited nodes, 0 failed targets. Heatmap showed:

  - Clear Voronoi-cell rays radiating from the origin (single-nearest-
    visited-node assignment changing abruptly between adjacent walk
    paths — discontinuities visible)
  - Spurious "poles" along the top/bottom edges at e.g. `(0, ±2.5)`
    (angle 90°, squarely in the analytic pole-free sector `|arg x| >
    36°` — so cannot be real poles)

The 8-target perimeter walk was too sparse: walks from the seed at
`-3 + 0im` to corners are √45 ≈ 6.7 units (vs √13 ≈ 3.6 for [-2, 2]²),
leaving large Voronoi cells where the Padé extrapolation is past its
convergence disc.

### Iteration 2 — 11×11 = 121 dense target lattice

Replaced the 8-perimeter targets with a 11×11 dense lattice over
`[-3, 3]²` at spacing 0.6.

**Result**: kernel 60.4 s, 1381 visited nodes (via `skip_covered = true`
+ adaptive stepping), 0 failed. Heatmap:

  - Voronoi-cell rays eliminated (the dense lattice fills the central
    region with overlapping walk paths)
  - Ghost poles at top/bottom edges PERSISTED — at slightly different
    positions and now with a regular spacing of ~1.0 along the edges

Diagnosis: even with denser targets, `extrapolate = true` evaluates the
nearest Padé past its convergence disc. Each Padé's `Q(z)` has zeros at
semi-random points in its analyticity region; these zeros do NOT
represent real solution poles, and when a fine-grid cell's nearest
visited-node has a `Q`-zero nearby, the cell renders as a spurious pole.
The parent surface kernel's VC-4 cross-node validation is what filters
these out; the zoom kernel skips VC-4 for simplicity.

### Iteration 3 — `extrapolate = false` + BFS nearest-NN fill

Changed `extrapolate = true → false` in `vector_path_network_solve` (the
zoom kernel) so Stage-2 cells past the nearest Padé's convergence disc
become `NaN`. Added a `_fill_zoom_nan!` BFS helper to both rendering
scripts (same idiom as the R=8 figures' `_fill_inner_disc_nan!` but with
no disc gate — the zoom is the full square).

**Result**: 49.2 % of cells `NaN` before fill, filled in 47 BFS
iterations to 100 % finite. Heatmap:

  - Ghost poles GONE (no past-disc Padé output anywhere)
  - REPLACED by visible polygonal BFS-fill regions at the perimeter —
    the BFS fill smoothed past-disc cells using values from the in-disc
    cells (which near pole regions are LEGITIMATELY high). The
    BFS-fill polygons against the dark smooth-sector background read
    as "blocky" pole-like features

The net visual quality dropped: ghost-pole-as-points became ghost-pole-
as-polygons. Iteration 3 wasn't a clean win.

### Iteration 4 — 21×21 = 441 targets

Doubled the target lattice density (spacing 0.6 → 0.3) while keeping
`extrapolate = false` + BFS fill.

**Result**: kernel 145.4 s, 3425 visited nodes (vs 1381), 0 failed.
NaN fraction dropped 49.2 % → 34.4 %; BFS iterations dropped 47 → 24.
Heatmap:

  - Polygonal BFS-fill regions thinned to discrete pole-like dots —
    visually MUCH cleaner, the pole field reads as a regular lattice
  - Left half (smooth sector toward `-3 + 0im`) cleanly dark
  - Right half (toward the wedge) shows a roughly regular pole lattice
  - Small textured patch around `(2.5, -1.5)` (insufficient coverage
    in that bottom-right corner)

This was a significant visual improvement. User-acknowledged.

### Iteration 5 — `h = 0.05` (walk step halved)

Halved the path-network walk step `h = 0.1 → 0.05` to tighten each
Padé's convergence disc (cleaner in-disc values, even if total covered
area drops).

**Result**: kernel 151.0 s, 3488 visited nodes — **only 2 % more than
h = 0.1**. The adaptive walk (`adaptive = true`, `:max_q_root`) was
already taking sub-`0.1` steps near pole regions; halving `h_max`
ceiling barely binds. NaN fraction 34.4 % → 32.9 %, BFS iterations
unchanged at 24. Heatmap was visually comparable to iteration 4 —
similar lattice, similar perimeter artefacts.

**Verdict**: diminishing returns. The structural limit isn't `h` or
target density; it's the single-nearest-node Stage-2 fill (one Padé per
cell) lacking cross-validation.

## What ships

The final configuration committed to `origin/main`:

  - `ZOOM_XY_LIM = 3.0` (the user-requested zoom over `[-3, 3]²`)
  - `ZOOM_GRID_N = 1001` (the user-requested resolution, `dx = 0.006`)
  - `target_lattice_n = 21` (441 dense lattice targets, spacing 0.3)
  - `ZOOM_PN_H = 0.05` (walk step ceiling, mostly binding away from
    poles where the adaptive controller doesn't shrink it)
  - `ZOOM_PN_ORDER = 48` (matches R=8 kernel)
  - `extrapolate = false` in the Stage-2 fill (cells past the nearest
    Padé's convergence disc become NaN)
  - `_fill_zoom_nan!` BFS nearest-non-NaN fill in BOTH zoom rendering
    scripts (papers over the 32.9 % past-disc NaN cells, 24 BFS
    iterations to 100 % finite)

Kernel runtime: 151 s (cache MISS); ~22 s (cache HIT, JLD2 deserialise +
Makie startup). Each rendering script ~30 s on cache HIT.

The two PNG files (`kkg_pi2_abs_heatmap_zoom.png` and
`kkg_pi2_abs_phase_surface_zoom.png`) are committed; the cache JLD2
stays gitignored.

## What's still imperfect (honest limitations)

  - **The bright spots at the perimeter** — most are likely real
    inner-wedge poles at `|x| ∈ [2, 3]`, but some may be BFS-fill-
    propagated values from cells whose Padé happened to have a
    `Q`-zero near the cell. Without K-nearest cross-validation we
    can't tell them apart.
  - **A small textured patch near `(2.5, -1.5)`** persists faintly in
    the heatmap — insufficient walk-path coverage in that bottom-right
    corner where the walks from seed travel furthest.
  - **The 3D phase surface** shows a forest of pole spikes; clusters at
    the perimeter blend "real near-wedge" and "artefact near-disc-
    boundary" features.

For the [-2, 2]² zoom (worklog 062) these limitations didn't matter
because the zoom region was entirely inside the analytic smooth sector
(no real poles) — the walk reached every cell with valid Padé output,
no BFS fill was needed, and the figure was unambiguous. The [-3, 3]²
zoom necessarily crosses the wedge inner boundary at `|x| = 2` and
includes pole-rich territory, which exposes the Stage-2 fill's
limitations.

## Pickup for the next agent — the principled fix

The structural limitation is **single-nearest-node Stage-2 fill**. The
fix is to evaluate the K nearest visited nodes' Padé approximants at
each fine-grid cell and average (or median) them. Spurious `Q`-zeros
from individual nodes get washed out; real poles (which all neighbours
agree on) survive.

This is `~1 day` of `src/` work in `src/VectorPathNetwork.jl`:

  - `_stage2_fill` (around line 815 of `vector_path_network_solve`)
    currently calls `_nearest_visited` for each grid cell.
  - Replace with `_k_nearest_visited(K)`; weight each contribution by
    `1 / (distance + ε)`.
  - Add a `k_nearest::Integer = 1` kwarg to `vector_path_network_solve`
    (default `1` preserves byte-identical existing behaviour).
  - Mutation-prove: with `k_nearest = 1` every existing test stays
    GREEN. With `k_nearest = 4` the zoom kernel's perimeter artefacts
    should soften substantially.

This is a P3-priority follow-up bead (a new bead `padetaylor-(k-nearest-
stage2)`); the current zoom figures are GREEN and ship in their honest-
limits state.

## Other follow-ups (carried from worklog 062)

  - Harmonise `_cache_is_fresh` across all four figure scripts (the
    dual-fill script's check is over-strict — includes the render-only
    script in the cache-invalidation set).
  - Wire `_load_or_compute_kernel` into `figures/test_kkg_pi2_surface.jl`
    so the VC suite shares the cache.
  - Optional: a `log10|f|` variant of the heatmaps for FW-paper-style
    figures (the smooth-sector detail would be more visible at the
    expense of the pole-saturation drama).

## Files touched this iteration arc

  - **`figures/_kkg_pi2_zoom_helpers.jl`** — `ZOOM_XY_LIM 2.0 → 3.0`,
    `ZOOM_GRID_N 801 → 1001`, `extrapolate = true → false`, 8
    perimeter targets → 21×21 dense lattice, `ZOOM_PN_H = 0.1 → 0.05`.
    Updated literate docstrings explaining each change.
  - **`figures/kkg_pi2_abs_heatmap_zoom.jl`** — added `_fill_zoom_nan!`
    BFS helper; updated docstring to describe the new fill mechanism;
    updated assertion to allow pre-fill NaN cells.
  - **`figures/kkg_pi2_abs_phase_surface_zoom.jl`** — same additions.
  - **`figures/output/kkg_pi2_abs_heatmap_zoom.png`** — replaced.
  - **`figures/output/kkg_pi2_abs_phase_surface_zoom.png`** — replaced.
  - **`.gitignore`** — added `figures/output/kkg_pi2_zoom_cache.jld2`.

## Commits

  - `d350ad5` tf9.5: zoom figures over [-2,2]² @ 801×801 (user-approved)
  - (this iteration arc — `[-3, 3]²` work — in commits to follow)

## References

  - worklog 062 — the parent figures' design
  - `figures/_kkg_pi2_zoom_helpers.jl` — the zoom kernel with full
    docstring rationale for each design choice
  - `src/VectorPathNetwork.jl:815` — `_stage2_fill` (the call site
    that would change for K-nearest cross-validation)
  - `figures/_kkg_pi2_vc45.jl:330` — `vc4_validate` (the cross-node
    validation idiom in the surface kernel that the zoom kernel
    deliberately skips)
  - beads `padetaylor-tf9.5` (closed) — the [-2, 2]² zoom
  - bead `padetaylor-(zoom-iteration-arc)` (this work)
