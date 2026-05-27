# figures/kkg_pi2_abs_phase_surface.jl
#
# A SINGLE-PANEL 3D MODULUS SURFACE, PHASE-COLOURED, of the
# P_I^(2) tritronquée `V₀(x, 0)` over the disc |x| ≤ SURF_XY_LIM (bead
# `padetaylor-tf9.4`).  Companion to the 2D `|V₀|` viridis heatmap in
# `kkg_pi2_abs_heatmap.jl`; both are *dishonest* renderings — they paint
# the FW-extrapolated kernel matrices (`res.Re_u_extrap`,
# `res.Im_u_extrap`) over the full disc with NO certified-vs-extrapolated
# provenance gate.  The B1-honest dual-fill 2×2 reference figure
# (`kkg_pi2_tritronquee_surface.{jl,png}`) is the provenance-marked
# headline; this figure is the standard *domain-colouring* view a
# complex-analysis textbook reader expects: `z = |f|` for the surface
# height, `colour = arg f` cyclic for the phase.
#
# ## Why this figure exists
#
# A single 3D `surface!` showing `|f|` AS HEIGHT and `arg f` AS COLOUR
# carries strictly more information than two separate Re/Im plots: the
# height encodes magnitude, the cyclic colour encodes phase, and pole
# structure reads as bright spikes around which the phase wraps a full
# `2π`.  This is the convention used by every modern complex-function
# explorer (Wegert's "Visual Complex Functions", `matplotlib`'s
# `colorize_complex`, `Mathematica`'s `ComplexPlot3D`).
#
# ## The "dishonest fill" — what gets papered over
#
# Same as the 2D heatmap companion: `res.Re_u_extrap` / `res.Im_u_extrap`
# are the FW-2011-style filled matrices (every Padé evaluated at every
# cell, no validity gate); the only NaN cells inside the disc are the
# `±1°` Stokes strips (`SURF_STITCH_MASK_DEG = 1.0` in
# `_kkg_pi2_surface_helpers.jl`).  `_fill_inner_disc_nan!` fills those
# strips by nearest-non-NaN-neighbour averaging so the surface reads as
# continuous — the user's "no empty spots" requirement.  Outside the
# disc the cells stay NaN; Makie's `surface!` automatically excludes
# NaN cells from the mesh, giving a clean disc-shaped surface boundary.
#
# ## Rendering choices
#
#   * **Height = `min(|f|, 15.0)`.**  `SURF_CLAMP = 15.0` matches the
#     dual-fill figure.  In the wedge `|f|` spikes to `O(10³)` near
#     poles; un-clamped, the 3D box would have a vertical aspect that
#     compresses the smooth sector to a flat plane.  Clamped, the pole
#     peaks read as a saturated plateau and the smooth sector retains
#     real relief.
#
#   * **Colour = `arg f`, cyclic colormap.**  We map `arg f ∈ [-π, π]`
#     to `[0, 1]` via `(arg f + π) / 2π` and feed it to a cyclic Makie
#     colormap.  Cyclic colormaps have no discontinuity at the ±π
#     wraparound — the phase is genuinely circular, and a non-cyclic
#     map (e.g. viridis) would draw a spurious seam along the negative
#     real axis.  Default cyclic choice: `:cyclic_mygbm_30_95_c78_n256`
#     (CET-C2 "mygbm" — magenta/yellow/green/blue/magenta) which has
#     equal perceptual lightness around the cycle so pole-encircling
#     phase wraps read as a continuous hue rotation, not a brightness
#     pulse.
#
#   * **Axis3 with `aspect = (1, 1, 0.5)`.**  The 3D box is wider in
#     Re/Im (the data plane) than tall in `|f|` (the clamped height).
#     `viewmode = :fit` lets Makie fit the box to the panel.
#
#   * **Phase colourbar with `arg` ticks at `±π, ±π/2, 0`.**  The
#     colourbar's `limits = (-π, π)` exposes the underlying phase scale
#     (not the [0,1] normalisation we feed `surface!`), so the reader
#     reads `arg V₀` directly.
#
# References:
#   * `references/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41.pdf`
#     — KKG 2015 §7 (the `V₀` solution).
#   * Wegert, E. "Visual Complex Functions" (2012) — the canonical
#     reference on phase-coloured complex surfaces.
#   * `figures/_kkg_pi2_surface_helpers.jl` — the kernel.
#   * `figures/kkg_pi2_tritronquee_surface.jl` — the B1-honest
#     dual-provenance dual-fill 2×2 reference figure.

using CairoMakie
using Printf
using JLD2

module KKGSurfaceKernel
include(joinpath(@__DIR__, "_kkg_pi2_surface_helpers.jl"))
end
using .KKGSurfaceKernel: kkg_pi2_surface, SURF_GRID_N, SURF_XY_LIM,
                         SURF_PN_ORDER

const OUTPNG        = joinpath(@__DIR__, "output", "kkg_pi2_abs_phase_surface.png")
const CACHE_PATH    = joinpath(@__DIR__, "output", "kkg_pi2_kernel_cache.jld2")
const CACHE_VERSION = "v1"
const SURF_CLAMP    = 15.0
const PHASE_CMAP    = :cyclic_mygbm_30_95_c78_n256

# --- Cache-aware kernel load ------------------------------------------
# Identical contract to `kkg_pi2_abs_heatmap.jl` — same JLD2 cache file,
# same header keys.  Render-only: the freshness gate excludes
# `basename(@__FILE__)` (no kernel logic here, so a new render script
# alone must not invalidate the cache).
function _kernel_cache_header()
    (version    = CACHE_VERSION,
     grid_n     = SURF_GRID_N,
     xy_lim     = SURF_XY_LIM,
     pn_order   = SURF_PN_ORDER,
     r_max      = KKGSurfaceKernel.SURF_R_MAX,
     pn_h       = KKGSurfaceKernel.SURF_PN_H)
end

function _cache_is_fresh(path)
    isfile(path) || return false
    haskey(ENV, "KKG_FORCE_RECOMPUTE") && return false
    cache_mtime = mtime(path)
    helper_dir  = @__DIR__
    for f in readdir(helper_dir; join = true)
        endswith(f, ".jl") || continue
        name = basename(f)
        if startswith(name, "_kkg_pi2")
            mtime(f) > cache_mtime && return false
        end
    end
    return true
end

function _load_or_compute_kernel()
    if _cache_is_fresh(CACHE_PATH)
        try
            cached_header, cached_res = JLD2.load(CACHE_PATH, "header", "res")
            if cached_header == _kernel_cache_header()
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
    res = kkg_pi2_surface()
    try
        mkpath(dirname(CACHE_PATH))
        JLD2.save(CACHE_PATH, "header", _kernel_cache_header(), "res", res)
        @printf("  cache: wrote %s\n", CACHE_PATH); flush(stdout)
    catch err
        @printf("  cache: write failed (%s) — continuing without persistence\n",
                err); flush(stdout)
    end
    return res
end

# --- Stokes-strip NaN fill --------------------------------------------
# Same helper as in `kkg_pi2_abs_heatmap.jl`.  See that file's docstring
# for the "dishonest visual fill" rationale.
function _fill_inner_disc_nan!(M::Matrix{Float64}, xs, ys, R::Real)
    n = size(M, 1)
    R2 = R * R
    iters = 0
    while true
        changed = false
        snap = copy(M)
        @inbounds for j in 1:n, i in 1:n
            isnan(M[i, j]) || continue
            xs[i]^2 + ys[j]^2 ≤ R2 || continue
            acc = 0.0; k = 0
            for dj in -1:1, di in -1:1
                (di == 0 && dj == 0) && continue
                ii = i + di; jj = j + dj
                (1 ≤ ii ≤ n && 1 ≤ jj ≤ n) || continue
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

# --- Main pipeline ----------------------------------------------------
@printf("phase: kkg_pi2_abs_phase_surface — 3D |V₀| height, phase-coloured\n"); flush(stdout)
@printf("  grid: %d×%d over [-%.0f,%.0f]² (R=%.0f)\n",
        SURF_GRID_N, SURF_GRID_N, SURF_XY_LIM, SURF_XY_LIM, SURF_XY_LIM); flush(stdout)

t0 = time()
@printf("phase: loading kernel (cache or compute)\n"); flush(stdout)
res = _load_or_compute_kernel()
@printf("phase: kernel ready in %.1f s — %s\n", time() - t0, res.message); flush(stdout)

xs, ys = res.xs, res.ys
Re_ext = copy(res.Re_u_extrap)
Im_ext = copy(res.Im_u_extrap)

@printf("phase: filling Stokes-strip NaN cells (Re_u_extrap)\n"); flush(stdout)
n_re_iters = _fill_inner_disc_nan!(Re_ext, xs, ys, SURF_XY_LIM)
@printf("  Re fill: %d BFS iterations\n", n_re_iters); flush(stdout)
@printf("phase: filling Stokes-strip NaN cells (Im_u_extrap)\n"); flush(stdout)
n_im_iters = _fill_inner_disc_nan!(Im_ext, xs, ys, SURF_XY_LIM)
@printf("  Im fill: %d BFS iterations\n", n_im_iters); flush(stdout)

# |f| height, clamped at SURF_CLAMP; arg f → [0,1] for the cyclic
# colormap input.  Both matrices share their NaN mask (the outside-disc
# cells); Makie's `surface!` skips NaN-z mesh nodes so the surface is
# automatically disc-clipped.
abs_clamped = [isnan(Re_ext[i, j]) || isnan(Im_ext[i, j]) ? NaN :
               min(sqrt(Re_ext[i, j]^2 + Im_ext[i, j]^2), SURF_CLAMP)
               for i in 1:size(Re_ext, 1), j in 1:size(Re_ext, 2)]
phase01     = [isnan(Re_ext[i, j]) || isnan(Im_ext[i, j]) ? NaN :
               (atan(Im_ext[i, j], Re_ext[i, j]) + π) / (2π)
               for i in 1:size(Re_ext, 1), j in 1:size(Re_ext, 2)]
@printf("  |V₀| clamped finite cells: %d / %d\n",
        count(isfinite, abs_clamped), length(abs_clamped)); flush(stdout)

# --- Render -----------------------------------------------------------
@printf("phase: building Makie figure\n"); flush(stdout)
t_render = time()
fig = Figure(size = (900, 780))
ax = Axis3(fig[1, 1];
           title = "|V₀|, phase-coloured  —  P_I⁽²⁾ tritronquée (FW-extrapolated)",
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

# Render-sanity check (CLAUDE.md Rule 1 — fail loud).
let bytes = filesize(OUTPNG)
    @assert isfile(OUTPNG) "render produced no PNG at $OUTPNG"
    @assert bytes > 10_000 "render PNG suspiciously small ($bytes bytes)"
    @printf("phase: render-sanity OK — PNG %d bytes\n", bytes); flush(stdout)
end
