# figures/kkg_pi2_abs_heatmap_zoom.jl
#
# A SINGLE-PANEL 2D HEATMAP of `|V₀(x, 0)|` for the P_I^(2) tritronquée
# over the **zoom window** `[-3, 3]²` at 1001 × 1001 resolution
# (`dx = 0.006`, ~6.7 × finer than the R = 8 cache).  Companion to
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
# integrates inward to every grid point of `[-3, 3]²`.  The seed is
# **on** the left edge of the zoom square (`Re = -3`), the bulk of the
# zoom is in the smooth pole-free sector (`|arg x| > 36°` covers most
# of `[-3, 3]²`), and the Stokes lines at `±36°` cross the zoom
# interior for `|x| ∈ [2, 3]` — the [-3, 3]² zoom shows wedge-sector
# poles in the band `|x| ∈ [2, 3], |arg x| < 36°`.  The single-walk
# kernel renders these poles as finite-amplitude spikes rather than
# separate wedge-walk artefacts.
#
# ## Stage-2 `extrapolate = false` + BFS fill
#
# At [-3, 3]² the seed-anchored walks travel up to √45 ≈ 6.7 units —
# well beyond the per-node Padé `h = 0.1` convergence disc.  Raw
# `extrapolate = true` Padé output past the disc produced SPURIOUS
# Q-roots (ghost poles) in the analytic smooth sector, user-flagged
# at the top/bottom edges (e.g. `(0, ±2.5)` at angle 90°, way outside
# the wedge `|arg x| < 36°` where the tritronquée has any real poles).
# The zoom kernel now runs Stage-2 with `extrapolate = false`: cells
# past the nearest visited Padé's disc become `NaN`.  A nearest-non-NaN
# BFS fill (`_fill_zoom_nan!` below — same Moore-neighbourhood idiom
# as the R = 8 `_fill_inner_disc_nan!` but with no disc gate, since
# the zoom IS the full square) papers those cells over with
# interpolated values from valid neighbours.  Result: the in-disc
# honest Padé evaluations stay intact, the past-disc cells are
# smoothed via interpolation, and the figure shows no ghost poles.
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

# Fill NaN cells from non-NaN neighbours via a Moore-3x3 BFS, iterated
# until stable.  Read-from-snapshot, write-to-live: one BFS layer per
# pass, no within-sweep cascading.  Same idiom as the R = 8 figure's
# `_fill_inner_disc_nan!` but with NO disc gate — the zoom is the
# full [-ZOOM_XY_LIM, ZOOM_XY_LIM]² square, no outside-disc cells.
function _fill_zoom_nan!(M::Matrix{Float64})
    n_i, n_j = size(M)
    iters = 0
    while true
        changed = false
        snap = copy(M)
        @inbounds for j in 1:n_j, i in 1:n_i
            isnan(M[i, j]) || continue
            acc = 0.0; k = 0
            for dj in -1:1, di in -1:1
                (di == 0 && dj == 0) && continue
                ii = i + di; jj = j + dj
                (1 ≤ ii ≤ n_i && 1 ≤ jj ≤ n_j) || continue
                v = snap[ii, jj]
                isnan(v) && continue
                acc += v; k += 1
            end
            if k > 0
                M[i, j] = acc / k
                changed = true
            end
        end
        iters += 1
        changed || break
    end
    return iters
end

Re_u = copy(res.Re_u)
Im_u = copy(res.Im_u)
n_nan_before = count(isnan, Re_u)
@printf("phase: BFS-filling past-disc NaN cells (Re_u then Im_u)\n"); flush(stdout)
@printf("  NaN cells in Re_u before fill: %d / %d (%.1f%%)\n",
        n_nan_before, length(Re_u),
        100.0 * n_nan_before / length(Re_u)); flush(stdout)
n_re_iters = _fill_zoom_nan!(Re_u)
n_im_iters = _fill_zoom_nan!(Im_u)
@printf("  Re_u fill: %d BFS iterations; Im_u fill: %d BFS iterations\n",
        n_re_iters, n_im_iters); flush(stdout)

abs_u = sqrt.(Re_u .^ 2 .+ Im_u .^ 2)
let n_finite = count(isfinite, abs_u)
    @printf("  |V₀| finite cells after fill: %d / %d\n",
            n_finite, length(abs_u)); flush(stdout)
    @assert n_finite == length(abs_u) (
        "zoom |V₀| has $(length(abs_u) - n_finite) NaN cells after BFS " *
        "fill — every cell should be filled.  Investigate: a NaN cell " *
        "with no non-NaN neighbours can only happen if a whole row/" *
        "column is NaN.")
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
