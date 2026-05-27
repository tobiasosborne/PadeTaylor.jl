# figures/kkg_pi2_abs_phase_surface_zoom.jl
#
# A SINGLE-PANEL 3D MODULUS SURFACE, PHASE-COLOURED, of `V₀(x, 0)`
# for the P_I^(2) tritronquée over the **zoom window** `[-2, 2]²` at
# 801 × 801 resolution.  Companion to `kkg_pi2_abs_phase_surface.jl`
# (the R = 8 version with the same `min(|f|, 15)` height, cyclic
# `:cyclic_mygbm_30_95_c78_n256` phase colormap, same `arg`-ticks
# colourbar).  The zoom shows the smooth-sector core of the surface
# with no wedge-walk artefacts — see the heatmap companion's docstring
# for the kernel-design rationale.
#
# ## Rendering choices — same as the R = 8 version
#
#   * **Height `z = min(|V₀|, 15.0)`.**  Matches `SURF_CLAMP = 15` from
#     the R = 8 figures.  In `|x| ≤ 2` (squarely in the smooth sector)
#     `|V₀|` is bounded by `O(10)` — the clamp barely engages near the
#     negative-real ridge floor and never engages in the calm bulk of
#     the zoom; the clamp is kept for consistency with the headline
#     figures rather than to tame any zoom-specific spike.
#
#   * **Colour `(arg V₀ + π) / 2π` through `:cyclic_mygbm_30_95_c78_n256`.**
#     The cyclic colormap renders the smooth-sector phase as a slow
#     hue rotation around the negative-real ridge (`arg V₀ ≈ 0` →
#     green-yellow at the ridge floor) with no spurious seam at the
#     `±π` wraparound.  Inside the zoom there are no poles to encircle
#     2π, so the phase varies smoothly — the colour gradient reads as
#     a continuous tinting of the surface rather than the discrete
#     pole-encircling pulses the R = 8 figure shows in its wedge.
#
#   * **`Axis3` with `aspect = (1, 1, 0.5)`, `viewmode = :fit`.**
#     Square in the data plane, half-height in the value axis — the
#     same proportions as the R = 8 figure.
#
#   * **Phase colourbar with `arg` ticks `±π, ±π/2, 0`.**  The
#     colourbar `limits = (-π, π)` exposes the underlying phase scale
#     (not the `[0, 1]` normalisation `surface!` consumes).
#
# References:
#   * `figures/_kkg_pi2_zoom_helpers.jl` — the zoom kernel.
#   * `figures/kkg_pi2_abs_phase_surface.jl` — the R = 8 companion.
#   * `figures/kkg_pi2_abs_heatmap_zoom.jl` — the 2D viridis companion;
#     shares the same `_load_or_compute_zoom_kernel` cache contract
#     (duplicated here so both render scripts are self-contained).

using CairoMakie
using Printf
using JLD2

module KKGZoomKernel
include(joinpath(@__DIR__, "_kkg_pi2_zoom_helpers.jl"))
end
using .KKGZoomKernel: kkg_pi2_zoom, ZOOM_XY_LIM, ZOOM_GRID_N, ZOOM_PN_ORDER

const OUTPNG        = joinpath(@__DIR__, "output", "kkg_pi2_abs_phase_surface_zoom.png")
const CACHE_PATH    = joinpath(@__DIR__, "output", "kkg_pi2_zoom_cache.jld2")
const CACHE_VERSION = "v1"
const SURF_CLAMP    = 15.0
const PHASE_CMAP    = :cyclic_mygbm_30_95_c78_n256

# --- Cache-aware kernel load ------------------------------------------
# Identical contract to `kkg_pi2_abs_heatmap_zoom.jl`.  Same JLD2 file,
# same header keys, same render-only freshness rule (only the zoom and
# surface helpers can invalidate the cache).
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
        cached_res = (xs = res.xs, ys = res.ys,
                      Re_u = res.Re_u, Im_u = res.Im_u,
                      n_visited = length(res.walk.visited_z),
                      n_failed  = length(res.walk.failed_targets),
                      message   = res.message)
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
@printf("phase: kkg_pi2_abs_phase_surface_zoom — 3D |V₀| height, phase-coloured\n"); flush(stdout)
@printf("  grid: %d×%d over [-%.1f, %.1f]², dx = %.4f\n",
        ZOOM_GRID_N, ZOOM_GRID_N, ZOOM_XY_LIM, ZOOM_XY_LIM,
        2 * ZOOM_XY_LIM / (ZOOM_GRID_N - 1)); flush(stdout)

t0 = time()
@printf("phase: loading kernel (cache or compute)\n"); flush(stdout)
res = _load_or_compute_zoom_kernel()
@printf("phase: kernel ready in %.1f s — %s\n", time() - t0, res.message); flush(stdout)

# Rule 1 (fail loud) — same invariant as the heatmap companion.
@assert res.n_failed == 0 (
    "zoom walk had $(res.n_failed) failed perimeter targets — investigate.")

xs, ys = res.xs, res.ys

# Build the height + phase fields.  The zoom kernel fills every cell,
# so the NaN-guard branch below is defensive (Rule 1: keep it; the
# branch is cheap and the assertion documents the invariant).
abs_clamped = [isnan(res.Re_u[i, j]) || isnan(res.Im_u[i, j]) ? NaN :
               min(sqrt(res.Re_u[i, j]^2 + res.Im_u[i, j]^2), SURF_CLAMP)
               for i in 1:size(res.Re_u, 1), j in 1:size(res.Re_u, 2)]
phase01     = [isnan(res.Re_u[i, j]) || isnan(res.Im_u[i, j]) ? NaN :
               (atan(res.Im_u[i, j], res.Re_u[i, j]) + π) / (2π)
               for i in 1:size(res.Re_u, 1), j in 1:size(res.Re_u, 2)]
let n_finite = count(isfinite, abs_clamped)
    @printf("  |V₀| clamped finite cells: %d / %d\n",
            n_finite, length(abs_clamped)); flush(stdout)
    @assert n_finite == length(abs_clamped) (
        "zoom surface has $(length(abs_clamped) - n_finite) NaN cells.")
end

# --- Render -----------------------------------------------------------
@printf("phase: building Makie figure\n"); flush(stdout)
t_render = time()
fig = Figure(size = (900, 780))
ax = Axis3(fig[1, 1];
           title = @sprintf("|V₀|, phase-coloured  —  P_I⁽²⁾ tritronquée zoom [-%.0f,%.0f]² @ %d×%d",
                            ZOOM_XY_LIM, ZOOM_XY_LIM, ZOOM_GRID_N, ZOOM_GRID_N),
           xlabel = "Re x", ylabel = "Im x", zlabel = "|V₀|",
           titlesize = 14,
           viewmode = :fit,
           aspect = (1.0, 1.0, 0.5))
surface!(ax, xs, ys, abs_clamped;
         color = phase01,
         colormap = PHASE_CMAP,
         colorrange = (0.0, 1.0),
         nan_color = RGBAf(1, 1, 1, 1))

Colorbar(fig[1, 2];
         colormap = PHASE_CMAP,
         limits = (-π, π),
         label = "arg(V₀)",
         ticks = ([-π, -π/2, 0.0, π/2, π],
                  ["-π", "-π/2", "0", "π/2", "π"]))

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("phase: wrote %s in %.1f s\n", OUTPNG, time() - t_render); flush(stdout)

# Render-sanity check (Rule 1 — fail loud).
let bytes = filesize(OUTPNG)
    @assert isfile(OUTPNG) "render produced no PNG at $OUTPNG"
    @assert bytes > 10_000 "render PNG suspiciously small ($bytes bytes)"
    @printf("phase: render-sanity OK — PNG %d bytes\n", bytes); flush(stdout)
end
