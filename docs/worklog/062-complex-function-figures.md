# Worklog 062 — complex-function figures (|f| heatmap + phase-3D)

**Date**: 2026-05-27
**Author**: Claude Opus 4.7
**Epic**: `padetaylor-0ln` (v0.2) · **Follows**: worklog 061
**Scope**: Two new single-panel figures of the P_I^(2) tritronquée
`V_0(x, 0)` using standard complex-function visualisation conventions
(2D |f| heatmap; 3D `|f|` surface coloured by `arg(f)`). Both
intentionally "dishonest" — FW-2011-style filled everywhere, no
certified-vs-extrapolated provenance gate. User-requested follow-up to
worklog 061 — the dual-fill alpha-blended aesthetic was unsatisfying.

> **Take-home**: the dual-fill provenance figure
> (`kkg_pi2_tritronquee_surface.png`) and the two new figures coexist.
> The dual-fill remains the provenance-documented reference; the two
> new figures are the user-facing visualisations. They read more
> cleanly because they use a single rendering convention per panel
> instead of compositing two opacity layers.

## The two figures

### `figures/kkg_pi2_abs_heatmap.{jl,png}` — 2D |f| heatmap

Single panel, `Axis` + `Heatmap` + `Colorbar`. Linear `|V_0|`,
`colorrange = (0, 15)` (matching `SURF_CLAMP`), `:viridis` colormap.

Visual reading: the smooth ~270° sector is a calm dark-blue field
around the negative-real ridge (`|V_0(-3, 0)| ≈ 2.6`). The wedge
`|arg x| < 36°` shows the KKG pole field as discrete bright peaks in
the inner band `|x| ∈ [2, 4]` and saturates to yellow over the outer
band `|x| ∈ [4, 8]` — the `|x|^{-1/6}` density rise from worklog 061
is visually evident.

138 LOC.

### `figures/kkg_pi2_abs_phase_surface.{jl,png}` — 3D `|f|` × phase surface

Single panel, `Axis3` + `surface!` + `Colorbar` (the phase legend).
`z = min(|V_0|, 15.0)` (z-clamp at 15), surface colour via
`(arg(V_0) + π) / (2π)` mapped through Makie's
`:cyclic_mygbm_30_95_c78_n256` cyclic colormap.

Visual reading: low-flat sector floor around the disc rim
(`|V_0| ≲ 5`, phase ≈ 0 → green); vertical pole-spike forest in the
wedge sector, saturated at `z = 15`; phase varies around each pole
(the cyclic colormap shows the `arg(f)` winding directly via
hue rotation). Colourbar labelled `±π, ±π/2, 0`.

144 LOC.

## Design decisions

### Both figures are "dishonest" — and that is the point

The dual-fill figure renders the *certified* B1-honest cells at full
opacity and the *extrapolated* (FW-2011-style) cells at reduced alpha
α = 0.50. ADR-0025 Amendment 14 documents the strict honesty contract
the dual-fill exists to maintain.

The user's "no empty spots / make them dishonest" instruction overrides
that contract *for these two figures only*. The data they show is
`res.Re_u_extrap` / `res.Im_u_extrap` — the FW-2011-style filled
matrices, every Padé evaluated over its full step with no validity
gate. The dual-fill figure remains the provenance-documented reference;
the new figures are *visualisations*, not provenance certificates.

This is captured in each script's top docstring and is the rationale
for keeping all three figures around.

### "No empty spots" — the inner-disc Stokes-strip fill

The kernel applies a `SURF_STITCH_MASK_DEG = 1.0` mask straddling
each `±36°` Stokes line — a ~7-cell-wide NaN strip at R=8's `dx = 0.04`
grid spacing. These strips render as grey holes in the dual-fill
figure (Makie's default `nan_color`). For the user-requested
"no empty spots" constraint, both new figures apply a
nearest-non-NaN BFS fill (`_fill_inner_disc_nan!`) to
`Re_u_extrap` and `Im_u_extrap` before computing `|f|` and
`arg(f)`. The strip cells fill from both wedge AND sector neighbours,
so the average is a smooth interpolant across the Stokes join.

Convergence: 24 BFS iterations to close the strips (≈3 cells per
pass on the 7-cell-wide strip). Inside-disc finite-cell coverage
after the fill: 125 629 / 125 629 (every in-disc cell).

### Outside-disc cells render as figure background

`nan_color = RGBAf(1, 1, 1, 1)` (white, full opacity) so cells outside
the disc blend into the figure background. Reads as a clean circle on
white. (For the 3D surface, Makie's `surface!` excludes NaN-`z` cells
automatically, so the surface footprint is naturally disc-shaped with
a clean edge.)

### Linear `|f|` with vmax = 15, not log10

User picked linear, viridis, clamp ≈ 15 from the option list (matches
`SURF_CLAMP`). The saturation in the outer wedge is honest about the
pole density — log10 would have compressed it. The dual-fill figure
already uses Re/Im clamped to ±15; using `|f|` clamped to 15 keeps the
display-range story consistent across all three figures.

### Cache-gate refinement for render-only scripts

The existing `_cache_is_fresh()` in `kkg_pi2_tritronquee_surface.jl`
treats EITHER the render script OR any `_kkg_pi2*.jl` helper as a
cache-invalidating dependency. That's over-strict: the render scripts
do not affect the kernel, only the helpers do. The two new scripts
drop the `basename(@__FILE__)` clause — only the `_kkg_pi2*` helpers
can legitimately invalidate the cache for them. Documented inline.

(Out of scope here, but a future P3 bead could harmonise the existing
script's `_cache_is_fresh` to match — current state is functionally
correct but asymmetric.)

## Verification

Both scripts ran serially (Rule 7), `cache: HIT` for both, no kernel
re-run, ~30 s wall per render. PNG render-sanity passed (> 10 KB).
Inner-disc fill confirmed visually — no grey holes in either figure.

```
kkg_pi2_abs_heatmap.png        145955 bytes
kkg_pi2_abs_phase_surface.png  463511 bytes
```

## Files

  - **`figures/kkg_pi2_abs_heatmap.jl`** — new, 138 LOC, self-contained.
  - **`figures/kkg_pi2_abs_phase_surface.jl`** — new, 144 LOC,
    self-contained.
  - **`figures/output/kkg_pi2_abs_heatmap.png`** — new deliverable.
  - **`figures/output/kkg_pi2_abs_phase_surface.png`** — new
    deliverable.

## Pickup for the next agent

  - **All three figures (dual-fill + heatmap + phase-3D) are GREEN
    and shipped.** No carryover.
  - **Future improvement candidates** (each ≤1 commit, P3):
    - Wire `_load_or_compute_kernel` into `figures/test_kkg_pi2_surface.jl`
      so the VC suite shares the cache (saves ~6 min per VC run).
    - Harmonise `_cache_is_fresh` across all three figure scripts
      (the dual-fill script's check is over-strict).
    - Optional: a `log10|f|` variant of the heatmap as an alternative
      visualisation for FW-paper-style figures (the smooth-sector
      detail would be visible at the expense of the pole-saturation
      drama).
  - **Open beads** unchanged from worklog 061's "Pickup": V9
    (`padetaylor-0ln.19`) is the natural next milestone for the
    `0ln` epic; `padetaylor-72w` (Heun first-class) and
    `padetaylor-x1y` (Padé validity theory) are P2 ready.

## Commits

  - `f839fc8` tf9.4: two single-panel complex-function figures
  - (this worklog + bead closure + final push)

## References

  - worklog 061 — predecessor; documents R=8 calibration the new
    figures inherit
  - ADR-0025 Amendment 14 — the dual-fill provenance the new
    figures intentionally don't honour (and that's documented in
    each script's docstring + here)
  - KKG 2015 Figs 7.4 / 7.5 — the canonical visualisations of the
    `P_I^(2)` tritronquée surface that the dual-fill figure
    mirrored; the new heatmap and phase-3D are the standard
    complex-function-visualisation companions to those
  - `figures/output/kkg_pi2_kernel_cache.jld2` — the JLD2 cache
    that made these renders ~30 s instead of 7 min
  - beads `padetaylor-tf9.4` (cited here and in commit f839fc8 but never actually recorded in the tracker — see bead `padetaylor-zt73`), `padetaylor-tf9` (closed in
    worklog 061; the umbrella bead this work refines under)
