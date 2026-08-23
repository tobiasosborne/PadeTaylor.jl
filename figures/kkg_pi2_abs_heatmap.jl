# figures/kkg_pi2_abs_heatmap.jl
#
# A SINGLE-PANEL 2D HEATMAP of the P_I^(2) tritronquée modulus |V₀(x,0)|
# over the disc |x| ≤ SURF_XY_LIM in the complex x-plane (bead
# `padetaylor-tf9.4` — an id that was never recorded in the tracker; the work is the closed parent `padetaylor-tf9`, see bead `padetaylor-zt73`).  Companion to the 3D modulus + phase surface in
# `kkg_pi2_abs_phase_surface.jl`; both are *dishonest* renderings — they
# present the FW-extrapolated kernel matrix (`res.Re_u_extrap`,
# `res.Im_u_extrap`) over the full disc with NO certified-vs-extrapolated
# provenance gate.  The B1-honest dual-fill 2×2 reference figure
# (`kkg_pi2_tritronquee_surface.{jl,png}`) is the provenance-marked
# headline; this figure is the FW-2011-style "fully filled" view a reader
# from outside the project intuits when they see "the complex-plane
# surface of the tritronquée".
#
# ## Why this figure exists
#
# KKG 2015 Figs 7.4/7.5 plot signed `Re V₀` and `Im V₀` — useful for
# verifying the sign-corrected IC (the negative-real-axis ridge must be
# positive, `u → +(6|x|)^{1/3}`) but un-intuitive at a glance.  The
# *standard* complex-function visualisation (e.g. Trefethen,
# `padeapprox.m`'s plotting harness, every Newton-fractal textbook
# figure) is `|f|` as a viridis heatmap: pole locations read as bright
# spots, the smooth analytic sector reads as a dark calm field, and the
# range of `|f|` is monotonic so a single colour-bar is unambiguous.
#
# ## The "dishonest fill" — what gets papered over
#
# The kernel returns two parallel matrices for each component:
#
#   * `res.Re_u` / `res.Im_u`  — the B1-honest surface: certified cells
#     (sector vote + wedge core whose Padé disc covers the point) and NaN
#     elsewhere.  ~29 % of the wedge is certified at order 48.
#   * `res.Re_u_extrap` / `res.Im_u_extrap` — the FW-2011-style surface:
#     every Padé evaluated at every cell, no validity gate.  Every
#     in-disc cell is finite EXCEPT the `±1°` Stokes strips
#     (`SURF_STITCH_MASK_DEG = 1.0` in `_kkg_pi2_surface_helpers.jl`),
#     which neither the sector fan reaches nor the wedge underlay
#     honestly covers.
#
# This figure paints the *extrapolated* matrices.  The Stokes strips are
# then filled by a small BFS-style nearest-neighbour averaging pass
# (`_fill_inner_disc_nan!`) so the disc reads as a continuous surface —
# the user's "no empty spots" requirement.  Outside the disc the cells
# stay NaN; Makie's `nan_color = :white` renders them flush with the
# figure background.
#
# ## Rendering choices
#
#   * **Colormap `:viridis`, range `(0, 15)`.**  Viridis is the standard
#     unsigned-field colormap (perceptually uniform, photocopier-safe).
#     The display clamp `SURF_CLAMP = 15.0` matches the dual-fill figure:
#     in the wedge `|V₀|` spikes to `O(10³)` near poles; without a clamp
#     the smooth sector (`|V₀| ≲ 3` along the negative-real ridge) would
#     wash to a flat single colour.  At 15, the smooth sector is fully
#     resolved and the clamped pole peaks read as a saturated plateau.
#
#   * **Square data aspect.**  The complex plane has equal Re/Im scales.
#
#   * **Colourbar label.**  `|V₀|  (display clamp ≤ 15)` — the clamp is
#     surfaced on the bar so a careful reader is not misled into reading
#     `|V₀| ≈ 15` as a true value (it's a colour saturation).
#
# References:
#   * `references/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41.pdf`
#     — KKG 2015 §7 Figs 7.4 (Re) / 7.5 (Im); the same `V₀` field, here
#     visualised as |·|.
#   * `references/markdown/FW2011_painleve_methodology_JCP230/`
#     — FW 2011 §5 published Padé pole-field figures fill without a
#     validity gate (FW-md:395-397); this script follows that convention.
#   * `figures/_kkg_pi2_surface_helpers.jl` — the kernel.
#   * `figures/kkg_pi2_tritronquee_surface.jl` — the B1-honest
#     dual-provenance dual-fill 2×2 reference figure (NOT replaced by
#     this one).

using CairoMakie
using Printf
using JLD2

# The kernel `include`s the Gridap-backed Laplace voter, whose
# `_kkg_pi2_gridap_helper.jl` uses an unqualified `Point`.  Pulled into
# the same `Main` as `using CairoMakie`, `Point` becomes ambiguous
# (CairoMakie/GeometryBasics vs Gridap.Fields both export it).  Isolate
# the kernel in its own module — identical pattern to the headline
# figure's `KKGSurfaceKernel`.
module KKGSurfaceKernel
include(joinpath(@__DIR__, "_kkg_pi2_surface_helpers.jl"))
end
using .KKGSurfaceKernel: kkg_pi2_surface, SURF_GRID_N, SURF_XY_LIM,
                         SURF_PN_ORDER

const OUTPNG       = joinpath(@__DIR__, "output", "kkg_pi2_abs_heatmap.png")
const CACHE_PATH   = joinpath(@__DIR__, "output", "kkg_pi2_kernel_cache.jld2")
const CACHE_VERSION = "v1"
const SURF_CLAMP   = 15.0

# --- Cache-aware kernel load ------------------------------------------
# Identical contract to `kkg_pi2_tritronquee_surface.jl`'s
# `_load_or_compute_kernel` — same JLD2 cache file, same `header` /
# `res` keys, same staleness rules.  Only the freshness check differs:
# THIS render script is render-only (no kernel logic), so it does NOT
# include `basename(@__FILE__)` in the mtime gate — only the `_kkg_pi2*`
# kernel helpers can invalidate the cache.  Otherwise a brand-new render
# script (newer mtime than the cache) would spuriously trigger a 3 h+
# recompute.
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
# Fill NaN cells INSIDE the disc by nearest non-NaN neighbour in the
# 3×3 Moore neighbourhood, iterated until stable.  This is a "dishonest"
# visual fill — it papers over the ±1° Stokes mask so the figure reads
# as the FW-2011-style filled surface the user wants.  Read-from-snapshot,
# write-to-live: each iteration is one BFS layer, no within-sweep
# cascading from cells written earlier in the same pass (which would
# bias the interior toward whichever scan direction we chose).
function _fill_inner_disc_nan!(M::Matrix{Float64}, xs, ys, R::Real)
    n = size(M, 1)
    R2 = R * R
    iters = 0
    while true
        changed = false
        snap = copy(M)
        @inbounds for j in 1:n, i in 1:n
            isnan(M[i, j]) || continue
            xs[i]^2 + ys[j]^2 ≤ R2 || continue   # outside disc — leave NaN
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
@printf("phase: kkg_pi2_abs_heatmap — single-panel |V₀| viridis heatmap\n"); flush(stdout)
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

# |f| = sqrt(Re² + Im²).  NaN cells outside the disc propagate naturally
# (NaN + NaN → NaN); inside-disc cells are all finite after the fill.
abs_u = sqrt.(Re_ext.^2 .+ Im_ext.^2)
let n_in_disc = count(i -> xs[((i-1) % size(abs_u,1))+1]^2 + ys[((i-1) ÷ size(abs_u,1))+1]^2 ≤ SURF_XY_LIM^2,
                     1:length(abs_u)),
    n_finite  = count(isfinite, abs_u)
    @printf("  |V₀| finite cells: %d / %d in-disc cells: %d\n",
            n_finite, length(abs_u), n_in_disc); flush(stdout)
end

# --- Render -----------------------------------------------------------
@printf("phase: building Makie figure\n"); flush(stdout)
t_render = time()
fig = Figure(size = (800, 720))
ax = Axis(fig[1, 1];
          title = "|V₀(x, 0)| — P_I⁽²⁾ tritronquée surface (FW-extrapolated)",
          xlabel = "Re x", ylabel = "Im x",
          aspect = DataAspect(),
          titlesize = 14,
          limits = (-SURF_XY_LIM, SURF_XY_LIM,
                    -SURF_XY_LIM, SURF_XY_LIM))
hm = heatmap!(ax, xs, ys, abs_u;
              colormap = :viridis,
              colorrange = (0.0, SURF_CLAMP),
              nan_color = RGBAf(1, 1, 1, 1))
Colorbar(fig[1, 2], hm;
         label = @sprintf("|V₀|  (display clamp ≤ %.0f)", SURF_CLAMP))

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("phase: wrote %s in %.1f s\n", OUTPNG, time() - t_render); flush(stdout)

# Render-sanity check (CLAUDE.md Rule 1 — fail loud).
let bytes = filesize(OUTPNG)
    @assert isfile(OUTPNG) "render produced no PNG at $OUTPNG"
    @assert bytes > 10_000 "render PNG suspiciously small ($bytes bytes)"
    @printf("phase: render-sanity OK — PNG %d bytes\n", bytes); flush(stdout)
end
