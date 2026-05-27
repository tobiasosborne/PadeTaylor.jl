# figures/kkg_pi2_abs_heatmap_zoom.jl
#
# A SINGLE-PANEL 2D HEATMAP of `|V₀(x, 0)|` for the P_I^(2) tritronquée
# over the **zoom window** `[-2, 2]²` at 801 × 801 resolution
# (`dx = 0.005`, 8 × finer than the R = 8 cache).  Companion to
# `kkg_pi2_abs_heatmap.jl` (the same field, viridis colormap, same
# `SURF_CLAMP = 15` saturation policy), but over the inner-disc smooth
# sector instead of the full `|x| ≤ 8` surface.
#
# ## Why the zoom is a separate kernel, not a re-render of the R = 8 cache
#
# The R = 8 figures stitch a sector-vote field onto a wedge-walk field
# across the `±36°` Stokes lines.  The visible artefacts the user
# reported — "false poles around `Im x = -2`, `Re x = 2.5`" — are in the
# inner-wedge / inner-sector join band at `|x| ≈ 2`–3.  The zoom kernel
# (`_kkg_pi2_zoom_helpers.jl`) bypasses that scaffolding entirely: a
# single Padé walk from the BVP-anchored seed at `z = -3 + 0im`
# integrates inward to every grid point of `[-2, 2]²`.  The seed is
# **outside** the zoom square (`Re = -3 < -2`), the bulk of the zoom is
# in the smooth pole-free sector (`|arg x| > 36°` covers all of
# `[-2, 2]²` except a small triangular tip at `(±2, 0)`), and the
# Stokes lines at `±36°` start at `|x| = 2` — they are at the boundary
# of the zoom, not inside it.  The result is one clean,
# wedge-artefact-free surface.
#
# ## The "dishonest fill" — same convention as the R = 8 figure
#
# The zoom kernel runs Stage-2 with `extrapolate = true`, so every grid
# cell receives a finite Padé evaluation regardless of the B1 true-
# radius gate.  This matches the worklog-062 convention the zoom
# inherits: the visible figure is FW-2011-style filled everywhere,
# and the headline-figure dual-fill provenance is preserved in the
# parallel `kkg_pi2_tritronquee_surface.png` figure.
#
# Unlike the R = 8 figure, the zoom needs NO `_fill_inner_disc_nan!`
# BFS step — the zoom kernel has no Stokes-strip mask (no `±36°`
# strips inside `|x| ≤ 2`), so every cell is finite on first pass.
# The render asserts this directly.
#
# ## Cache contract
#
# Cache file: `figures/output/kkg_pi2_zoom_cache.jld2` (separate from
# the R = 8 cache — different grid, different solver topology).
# Freshness: only the zoom or surface helper-file mtime invalidates
# the cache.  The render script itself does NOT invalidate (the
# worklog-062 refinement: render-only scripts have no kernel logic).
#
# References:
#   * `figures/_kkg_pi2_zoom_helpers.jl` — the zoom kernel.
#   * `figures/kkg_pi2_abs_heatmap.jl` — the R = 8 companion this
#     figure complements; same viridis colormap, same `SURF_CLAMP = 15`
#     saturation policy.
#   * `docs/worklog/062-complex-function-figures.md` — the parent
#     figures' design notes.

using CairoMakie
using Printf
using JLD2

# Isolate the kernel in its own module — the surface helpers
# transitively pull in Gridap, which exports an unqualified `Point`
# that clashes with CairoMakie's.  Same pattern as the R = 8 scripts.
module KKGZoomKernel
include(joinpath(@__DIR__, "_kkg_pi2_zoom_helpers.jl"))
end
using .KKGZoomKernel: kkg_pi2_zoom, ZOOM_XY_LIM, ZOOM_GRID_N, ZOOM_PN_ORDER

const OUTPNG        = joinpath(@__DIR__, "output", "kkg_pi2_abs_heatmap_zoom.png")
const CACHE_PATH    = joinpath(@__DIR__, "output", "kkg_pi2_zoom_cache.jld2")
const CACHE_VERSION = "v1"
const SURF_CLAMP    = 15.0

# --- Cache-aware kernel load ------------------------------------------
# Same JLD2 cache contract as `kkg_pi2_abs_heatmap.jl`: header captures
# every kernel-defining constant, mismatch triggers recompute.  Render-
# only freshness gate (no `basename(@__FILE__)` in the mtime check) —
# only `_kkg_pi2_zoom_helpers.jl` and `_kkg_pi2_surface_helpers.jl`
# (the latter is `include`d by the former for `surf_anchor_bvp` re-use)
# can legitimately invalidate the cache.
function _zoom_cache_header()
    (version  = CACHE_VERSION,
     grid_n   = ZOOM_GRID_N,
     xy_lim   = ZOOM_XY_LIM,
     pn_order = ZOOM_PN_ORDER,
     pn_h     = KKGZoomKernel.ZOOM_PN_H,
     pn_tol   = KKGZoomKernel.ZOOM_PN_TOL,
     t        = KKGZoomKernel.ZOOM_T)
end

function _zoom_cache_is_fresh(path)
    isfile(path) || return false
    haskey(ENV, "KKG_FORCE_RECOMPUTE") && return false
    cache_mtime = mtime(path)
    helper_dir  = @__DIR__
    for f in readdir(helper_dir; join = true)
        endswith(f, ".jl") || continue
        name = basename(f)
        if name == "_kkg_pi2_zoom_helpers.jl" ||
           name == "_kkg_pi2_surface_helpers.jl"
            mtime(f) > cache_mtime && return false
        end
    end
    return true
end

function _load_or_compute_zoom_kernel()
    if _zoom_cache_is_fresh(CACHE_PATH)
        try
            cached_header, cached_res = JLD2.load(CACHE_PATH, "header", "res")
            if cached_header == _zoom_cache_header()
                @printf("  cache: HIT (%s)\n", CACHE_PATH); flush(stdout)
                return cached_res
            end
            @printf("  cache: header mismatch, recomputing\n"); flush(stdout)
        catch err
            @printf("  cache: load failed (%s), recomputing\n", err); flush(stdout)
        end
    else
        @printf("  cache: MISS (%s); KKG_FORCE_RECOMPUTE=%s\n",
                CACHE_PATH, get(ENV, "KKG_FORCE_RECOMPUTE", "unset"));
        flush(stdout)
    end
    res = kkg_pi2_zoom()
    try
        mkpath(dirname(CACHE_PATH))
        # The `walk` field is bulky (every visited node's jets +
        # numerators + denominator), only useful for diagnostics.
        # Strip it from the cache so reloads are fast and the JLD2
        # file stays manageable.
        cached_res = (xs = res.xs, ys = res.ys,
                      Re_u = res.Re_u, Im_u = res.Im_u,
                      n_visited      = length(res.walk.visited_z),
                      n_failed       = length(res.walk.failed_targets),
                      message        = res.message)
        JLD2.save(CACHE_PATH, "header", _zoom_cache_header(), "res", cached_res)
        @printf("  cache: wrote %s\n", CACHE_PATH); flush(stdout)
        return cached_res
    catch err
        @printf("  cache: write failed (%s) — continuing without persistence\n",
                err); flush(stdout)
        return (xs = res.xs, ys = res.ys,
                Re_u = res.Re_u, Im_u = res.Im_u,
                n_visited = length(res.walk.visited_z),
                n_failed  = length(res.walk.failed_targets),
                message   = res.message)
    end
end

# --- Main pipeline ----------------------------------------------------
@printf("phase: kkg_pi2_abs_heatmap_zoom — |V₀| viridis heatmap over [-%.1f, %.1f]²\n",
        ZOOM_XY_LIM, ZOOM_XY_LIM); flush(stdout)
@printf("  grid: %d×%d, dx = %.4f\n",
        ZOOM_GRID_N, ZOOM_GRID_N, 2 * ZOOM_XY_LIM / (ZOOM_GRID_N - 1)); flush(stdout)

t0 = time()
@printf("phase: loading kernel (cache or compute)\n"); flush(stdout)
res = _load_or_compute_zoom_kernel()
@printf("phase: kernel ready in %.1f s — %s\n", time() - t0, res.message); flush(stdout)

# Rule 1 (fail loud): no failed perimeter targets allowed.  A non-zero
# count means the seed walk from `z = -3` could not reach a corner of
# the zoom square; that is a physical signal the kernel needs to
# investigate, not silently fill.
@assert res.n_failed == 0 (
    "zoom walk had $(res.n_failed) failed perimeter targets — the seed " *
    "walk from z = SURF_Z_SEED could not reach every corner of the " *
    "[-$(ZOOM_XY_LIM), $(ZOOM_XY_LIM)]² square.  Investigate (do NOT " *
    "relax on_target_failure).")

xs, ys = res.xs, res.ys

# `|V₀|` is finite everywhere — the zoom kernel has no Stokes-strip
# mask, no provenance NaN, no outside-disc clipping (the zoom *is* a
# square, not a disc).
abs_u = sqrt.(res.Re_u .^ 2 .+ res.Im_u .^ 2)
let n_finite = count(isfinite, abs_u)
    @printf("  |V₀| finite cells: %d / %d\n",
            n_finite, length(abs_u)); flush(stdout)
    @assert n_finite == length(abs_u) (
        "zoom |V₀| has $(length(abs_u) - n_finite) NaN cells — the zoom " *
        "kernel should fill every cell.  Investigate.")
end

# --- Render -----------------------------------------------------------
@printf("phase: building Makie figure\n"); flush(stdout)
t_render = time()
fig = Figure(size = (800, 720))
ax = Axis(fig[1, 1];
          title = @sprintf("|V₀(x, 0)|  —  P_I⁽²⁾ tritronquée zoom [-%.0f,%.0f]² @ %d×%d",
                           ZOOM_XY_LIM, ZOOM_XY_LIM, ZOOM_GRID_N, ZOOM_GRID_N),
          xlabel = "Re x", ylabel = "Im x",
          aspect = DataAspect(),
          titlesize = 14,
          limits = (-ZOOM_XY_LIM, ZOOM_XY_LIM,
                    -ZOOM_XY_LIM, ZOOM_XY_LIM))
hm = heatmap!(ax, xs, ys, abs_u;
              colormap = :viridis,
              colorrange = (0.0, SURF_CLAMP),
              nan_color = RGBAf(1, 1, 1, 1))
Colorbar(fig[1, 2], hm;
         label = @sprintf("|V₀|  (display clamp ≤ %.0f)", SURF_CLAMP))

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("phase: wrote %s in %.1f s\n", OUTPNG, time() - t_render); flush(stdout)

# Render-sanity check (Rule 1 — fail loud).
let bytes = filesize(OUTPNG)
    @assert isfile(OUTPNG) "render produced no PNG at $OUTPNG"
    @assert bytes > 10_000 "render PNG suspiciously small ($bytes bytes)"
    @printf("phase: render-sanity OK — PNG %d bytes\n", bytes); flush(stdout)
end
