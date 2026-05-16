# figures/ffw2017_fig_3.jl
#
# Reproduces Fasondini, Fornberg & Weideman 2017, **Figure 3** — the
# PVI phase-portrait reveal of the tronquée pole-free sector on the
# Riemann surface in BOTH the ζ-plane and the z-plane.
#
# Source: references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
#         FFW2017_painleve_riemann_surfaces_preprint.md:178-199.  The
#         caption verbatim (md:199):
#
# > Figure 3: The phase of the PVI solution shows the structure of the
# > Riemann surfaces in the ζ and z planes (recall that z = e ζ).  The
# > phase of the solution is depicted according to the color wheel,
# > taken from <http://dlmf.nist.gov/help/vrml/aboutcolor>.  The
# > pole-free 0th sheet in the z-plane only has a branch cut on
# > z ∈ (0, 1) whereas the other sheets have branch cuts on z ∈ (0, 1)
# > and the negative real z-axis.
#
# Same PVI problem as Figure 2 (FFW md:195 verbatim parameters + IC):
# `(α, β, γ, δ) = (4, -4, 8, -8)` with `u(10) = 0.429534600325223`,
# `u'(10) = -1.61713114374804e-3`.  Same A4 cross-mode + A5 sheet-aware
# Stage-2 walker as Figure 2's column 3.
#
# ## What this figure adds over Fig 2
#
# Where Fig 2 shows `|w(ζ)|` modulus heatmaps across three METHODS for
# the same solution (one column per method), Fig 3 shows the
# **phase-portrait + modulus** of the cross-mode (A4+A5) solution on
# TWO sheets in BOTH coordinate frames simultaneously:
#
#   * Top row — sheet [0] (principal): the pole-free 0th sheet.  The
#     z-plane branch cut is `z ∈ (0, 1)` only.
#   * Bottom row — alternate sheet (sheet bumped by the walker
#     crossing the ζ = 0 cut): branch cuts both on `z ∈ (0, 1)` AND
#     along the negative real z-axis.
#
# Three columns per row:
#
#   1. `|w(ζ)|` modulus heatmap on the ζ-plane (viridis, capped).
#   2. `arg w(ζ)` phase portrait on the ζ-plane (HSV color wheel).
#      The branch cut along `Re ζ < 0` (the negative-real ray emanating
#      from ζ = 0) is unmistakable on the bottom row: the phase jumps
#      across the cut.
#   3. `arg u(z)` phase portrait on the z-plane (HSV color wheel).
#      Computed via Cartesian-resample bilinear interpolation from the
#      ζ-grid `W` array — the same pattern as `figures/ffw2017_fig_4.jl`
#      (worklog 038 §"Column 2 z-plane Cartesian-resample bilinear-interp").
#      Sheet `s`: ζ = log|z| + i (arg z + 2π·s).
#
# ## ζ-window choice (v1 corner per CLAUDE.md Rule 9)
#
# FFW Fig 3's z-window is `|z| ∈ [1/100, 10]` (md:183), i.e.
# `Re ζ ∈ [log(1/100), log(10)] ≈ [-4.61, 2.30]`.  Our v1 ships a
# **reduced window** `Re ζ ∈ [-1.5, 2.5]` (`|z| ∈ [0.22, 12.2]`) —
# wider than Fig 2's `[0.05, 3.0]` but not blowing up walker wall-time
# on the low-Re ζ side, where the pole density grows rapidly close to
# the ζ = 0 branch point (FFW md:178 "high pole densities in the
# neighborhoods of the branch points").  The narrower window keeps the
# walker stable and the figure-script wall ≤ 30 s.  Extending to FFW's
# full `[-4.6, 2.3]` window is a deferred upgrade — same kind of
# Poisson-disk-node lift as bead `padetaylor-zwh`.
#
# ## rng_seed sensitivity (honest v1 corner)
#
# The walker's `shuffle(rng, targets)` step (FW 2011 line 156) creates
# rng-dependent visited trees.  For our wider Im-window window
# (`Im ζ ∈ (-π, 3π]` — same as Fig 2 col 3) + winding ring of 8 targets,
# `rng_seed = 0` (the package default) happens to produce a target
# order under which the walker never crosses the negative-real cut at
# `Re ζ < 0`: all 247 visited nodes stay on sheet [0].  A small sweep
# (probed during figure development) shows seeds 1, 2, 17, 42 all
# populate the alternate sheet; we use `rng_seed = 2` (526 visited
# nodes, sheet distribution `[0] => 477, [1] => 49`).  This is rng
# fragility we own honestly per CLAUDE.md Rule 9 — the principled fix
# is a Poisson-disk target placement that naturally circumambulates
# the branch point (same deferred infrastructure as bead
# `padetaylor-zwh`).
#
# ## Walker-targets vs render-lattice split (Fig 2 / Fig 6 pattern)
#
# Per `figures/ffw2017_fig_2.jl` header, `path_network_solve` walks the
# walker to EVERY entry in `grid`.  We supply ~150 SPARSE Stage-1
# targets (rectangular base + 8 winding-ring targets straddling
# `Re ζ < 0` to force the walker through the cut at least once),
# then evaluate Stage-2 manually on a dense render lattice via
# `eval_at_sheet(sol, z, sheet; extrapolate=true)`.  `extrapolate=true`
# (ADR-0015) fills the rendered panels without white gaps at the cost
# of Padé-extrapolation past `|t| = 1` for cells far from any visited
# node.

using PadeTaylor
using PadeTaylor.SheetTracker:   pVI_transformed_rhs
using PadeTaylor.CoordTransforms: pV_z_to_ζ, pV_ζ_to_z
using CairoMakie
using Printf

# ----------------------------------------------------------------------
# Parameters and ICs (FFW md:195 verbatim — same as Fig 2)
# ----------------------------------------------------------------------
const α, β, γ, δ = 4.0, -4.0, 8.0, -8.0
const Z0, U0, UP0 = 10.0+0im, 0.429534600325223+0im, -1.61713114374804e-3+0im

const ORDER       = 30
const ADAPT_TOL   = 1.0e-10
const K_CONS      = 1.0e-3
const MAX_RESC    = 50

# ζ-frame IC via the PV exponential transform (PVI reuses PV's exp
# coord map; FFW md:146).
const ζ0, w0, wp0 = pV_z_to_ζ(Z0, U0, UP0)

@printf("FFW 2017 Fig 3: (α, β, γ, δ) = (%.1f, %.1f, %.1f, %.1f)\n",
        α, β, γ, δ)
@printf("                IC (z₀, u, u') = (%s, %s, %s)\n", Z0, U0, UP0)
@printf("                  → (ζ₀, w, w') ≈ (%s, %s, %s)\n", ζ0, w0, wp0)

# ----------------------------------------------------------------------
# Window definitions
# ----------------------------------------------------------------------
# ζ-window: wider than Fig 2 (which uses Re ζ ∈ [0.05, 3.0]).  FFW's
# full Fig 3 window is `Re ζ ∈ [-4.6, 2.3]`; we ship the central band
# [-1.5, 2.5] as the v1 corner (see docstring).
const ζ_RE_LO, ζ_RE_HI   = -1.5, 2.5
const ζ_IM_LO, ζ_IM_HI   = -π + 0.05, 3π - 0.05

# Stage-2 render lattice density (rectangular ζ-grid; z-plane panels
# derived by Cartesian-resample bilinear interpolation).
const NX, NY = 80, 200

# z-plane render half-window matches |z| ≤ exp(Re ζ_HI) ≈ 12.2.
const Z_RMAX = exp(ζ_RE_HI)
const Z_HALF = 1.05 * Z_RMAX
const Z_NX, Z_NY = 280, 280

# ----------------------------------------------------------------------
# Build SPARSE rectangular Stage-1 target grid (NOT the render lattice).
# Same `0.35 × 0.7` spacing as Fig 2's cross-mode column.  Adds the
# 8 winding-ring targets that force the walker to cross the cut.
# ----------------------------------------------------------------------
function rect_targets(re_lo, re_hi, im_lo, im_hi, dx, dy)
    tgs = ComplexF64[]
    re = re_lo
    while re ≤ re_hi + 1e-12
        im_ = im_lo
        while im_ ≤ im_hi + 1e-12
            push!(tgs, complex(re, im_))
            im_ += dy
        end
        re += dx
    end
    return tgs
end

# Base targets: positive-Re region ONLY.  If we let the rectangular
# raster include `Re ∈ [-1.5, 0]`, the walker walks from positive Re
# DOWN through positive Im to reach those negative-Re targets without
# ever crossing `y = 0` at `Re < 0` — i.e. without crossing the cut
# (and visited_sheet stays at [0] throughout, just like dropping
# cross_branch=true).  Restricting base to Re > 0 forces the negative-
# Re cells to be reached only via the winding ring's CUT-CROSSING
# paths.
const targets_base = rect_targets(0.05, ζ_RE_HI, ζ_IM_LO, ζ_IM_HI,
                                   0.35, 0.7)
# Winding ring: pairs of upper-half and lower-half targets at Re < 0.
# The walker, having reached a (-1.0, +1.0i)-style target via upper-
# half travel, then has to reach (-1.0, -1.0i) — the min-|u| path
# generally cuts back through y = 0 at Re < 0, bumping the sheet
# counter.  Denser ring than Fig 2 because Fig 3's window extends to
# Re ζ = -1.5.
const winding_ring = ComplexF64[
    -0.5 + 0.5im, -0.5 - 0.5im,
    -1.0 + 0.5im, -1.0 - 0.5im,
    -1.0 + 1.0im, -1.0 - 1.0im,
    -0.5 + 1.0im, -0.5 - 1.0im,
]
const stage1_targets = vcat(targets_base, winding_ring)
@printf("Stage-1 targets: %d (incl. %d winding ring)\n",
        length(stage1_targets), length(winding_ring))

# ----------------------------------------------------------------------
# Build the ζ-frame PVI problem and walk it.
# ----------------------------------------------------------------------
const f_ζ = pVI_transformed_rhs(α, β, γ, δ)
const prob_ζ = PadeTaylorProblem(f_ζ, (w0, wp0),
                                  (ζ0, complex(ζ_RE_HI, ζ_IM_HI));
                                  order = ORDER)

println("\n== Stage 1: A4 cross-mode walk + A5 sheet bookkeeping ==")
t0 = time()
const RNG_SEED = 2
const sol = path_network_solve(prob_ζ, stage1_targets;
                                h = 0.3,
                                step_size_policy = :adaptive_ffw,
                                adaptive_tol = ADAPT_TOL,
                                k_conservative = K_CONS,
                                max_rescales = MAX_RESC,
                                max_steps_per_target = 4000,
                                branch_points     = (0.0 + 0.0im,),
                                branch_cut_angles = π,
                                cross_branch      = true,
                                rng_seed          = RNG_SEED)
@printf("  Stage 1 in %.2f s; %d visited tree nodes\n",
        time() - t0, length(sol.visited_z))

const sheet_counts = Dict{Vector{Int},Int}()
for s in sol.visited_sheet
    sheet_counts[s] = get(sheet_counts, s, 0) + 1
end
println("  Visited-sheet distribution: ", sheet_counts)

# ----------------------------------------------------------------------
# Pick the alternate sheet: the non-[0] sheet with the most visited
# nodes (same heuristic as Fig 2 §"Alternate-sheet rendering target").
# ----------------------------------------------------------------------
const alt_sheet_label = let
    nonzero = [(s, n) for (s, n) in sheet_counts if s != [0]]
    isempty(nonzero) ? Int[] :
        first(sort!(nonzero; by = x -> -x[2]))[1]
end
@printf("  Alternate-sheet rendering target: %s\n", alt_sheet_label)

# ----------------------------------------------------------------------
# Stage-2 render lattices.  ζ-plane: full-window rectangular grid;
# render two sheet panels by per-cell `eval_at_sheet`.
# ----------------------------------------------------------------------
const ζ_xs = collect(range(ζ_RE_LO, ζ_RE_HI; length = NX))
const ζ_ys = collect(range(ζ_IM_LO, ζ_IM_HI; length = NY))

W_sheet0 = Matrix{ComplexF64}(undef, NX, NY)
W_sheetalt = Matrix{ComplexF64}(undef, NX, NY)
t0_eval = time()
for j in 1:NY, i in 1:NX
    z = complex(ζ_xs[i], ζ_ys[j])
    u0_at, _ = eval_at_sheet(sol, z, [0]; extrapolate = true)
    W_sheet0[i, j] = u0_at
    if !isempty(alt_sheet_label)
        ualt, _ = eval_at_sheet(sol, z, alt_sheet_label;
                                extrapolate = true)
        W_sheetalt[i, j] = ualt
    else
        W_sheetalt[i, j] = complex(NaN, NaN)
    end
end
@printf("  ζ-plane Stage-2 eval in %.2f s\n", time() - t0_eval)

const cov_s0  = count(isfinite, abs.(W_sheet0))
const cov_alt = count(isfinite, abs.(W_sheetalt))
@printf("  ζ Stage-2 coverage: sheet0 %d/%d (%.1f%%), sheet%s %d/%d (%.1f%%)\n",
        cov_s0, length(W_sheet0), 100*cov_s0/length(W_sheet0),
        alt_sheet_label, cov_alt, length(W_sheetalt),
        100*cov_alt/length(W_sheetalt))

# ----------------------------------------------------------------------
# Bilinear interpolation of `w(ζ)` from the regular ζ-grid (same as
# `figures/ffw2017_fig_4.jl`'s `bilinear_w`).  Drives the z-plane
# panel rendering — no polar-scatter spoke artefacts.
# ----------------------------------------------------------------------
function bilinear_w(W::Matrix{ComplexF64}, ζ::Complex)
    rx, iy_ = real(ζ), imag(ζ)
    (rx < ζ_xs[1] || rx > ζ_xs[end] ||
     iy_ < ζ_ys[1] || iy_ > ζ_ys[end]) &&
        return complex(NaN, NaN)
    dx = ζ_xs[2] - ζ_xs[1]
    dy = ζ_ys[2] - ζ_ys[1]
    i = clamp(floor(Int, (rx - ζ_xs[1]) / dx) + 1, 1, NX - 1)
    j = clamp(floor(Int, (iy_ - ζ_ys[1]) / dy) + 1, 1, NY - 1)
    tx = (rx - ζ_xs[i]) / dx
    ty = (iy_ - ζ_ys[j]) / dy
    w00 = W[i, j]; w10 = W[i+1, j]
    w01 = W[i, j+1]; w11 = W[i+1, j+1]
    (isfinite(real(w00)) && isfinite(real(w10)) &&
     isfinite(real(w01)) && isfinite(real(w11))) ||
        return complex(NaN, NaN)
    return (1 - tx) * (1 - ty) * w00 +
           tx       * (1 - ty) * w10 +
           (1 - tx) * ty       * w01 +
           tx       * ty       * w11
end

# ----------------------------------------------------------------------
# z-plane render: for each sheet `s ∈ {0, 1}`, compute Cartesian z grid
# of |u(z)| + arg u(z).  PV transform (PVI inherits it): `u(z) = w(ζ)`
# with `ζ = log|z| + i (arg z + 2π · s)`.  Sheet 0 uses principal
# arg ∈ (-π, π]; sheet 1 lifts arg by 2π.
# ----------------------------------------------------------------------
const z_xs = collect(range(-Z_HALF, Z_HALF; length = Z_NX))
const z_ys = collect(range(-Z_HALF, Z_HALF; length = Z_NY))

# We render the alternate sheet's z-plane as the "sheet 1" face — even
# if the walker bumped to sheet [-1], the z-plane projection only sees
# `arg z + 2π · 1` vs `arg z + 0`; either bumps the cut convention the
# same way for the purpose of the phase-portrait reveal.
const alt_sign = (isempty(alt_sheet_label) || first(alt_sheet_label) == 0) ?
                  1 : sign(first(alt_sheet_label))

function render_zplane(W::Matrix{ComplexF64}, sheet_lift::Int)
    U_abs = fill(NaN, Z_NX, Z_NY)
    U_phi = fill(NaN, Z_NX, Z_NY)
    for jx in 1:Z_NY, ix in 1:Z_NX
        z = complex(z_xs[ix], z_ys[jx])
        absz = abs(z)
        absz < exp(ζ_RE_LO) - 0.005 && continue
        absz > exp(ζ_RE_HI) + 0.005 && continue
        ζcell = complex(log(absz), angle(z) + 2π * sheet_lift)
        # If lifted ζ falls outside our walked window in Im ζ, skip.
        (imag(ζcell) < ζ_IM_LO || imag(ζcell) > ζ_IM_HI) && continue
        w = bilinear_w(W, ζcell)
        isfinite(real(w)) || continue
        # PV / PVI: u(z) = w(ζ) (no rescale).
        U_abs[ix, jx] = abs(w)
        U_phi[ix, jx] = Float64(angle(w))
    end
    return U_abs, U_phi
end

t0_z = time()
const Uabs_z_s0, Uphi_z_s0   = render_zplane(W_sheet0,   0)
const Uabs_z_alt, Uphi_z_alt = render_zplane(W_sheetalt, alt_sign)
@printf("  z-plane render in %.2f s\n", time() - t0_z)

# ----------------------------------------------------------------------
# Render — 2 rows × 3 cols.
#   Col 1: |w(ζ)| modulus heatmap (viridis)
#   Col 2: arg w(ζ) phase portrait (HSV)
#   Col 3: arg u(z) phase portrait on z-plane (HSV)
# Top row: sheet [0].  Bottom row: alternate sheet.
# ----------------------------------------------------------------------
const U_CAP = 6.0
clip(x) = isfinite(x) ? min(x, U_CAP) : NaN

const OUTPNG = joinpath(@__DIR__, "output", "ffw2017_fig_3.png")
mkpath(dirname(OUTPNG))

fig = Figure(size = (1700, 1050))
Label(fig[0, 1:3],
      "FFW 2017 Fig 3 — PVI phase portraits on ζ + z Riemann surfaces, " *
      "(α,β,γ,δ) = $((α, β, γ, δ)), u(10) = $U0";
      fontsize = 15, padding = (0, 0, 4, 4))

# --- Top row: sheet [0] -----------------------------------------------
ax_abs0 = Axis(fig[1, 1];
               title = "(1a) |w(ζ)|, sheet [0]",
               xlabel = "Re ζ", ylabel = "Im ζ",
               aspect = DataAspect(),
               limits = (ζ_RE_LO, ζ_RE_HI, ζ_IM_LO, ζ_IM_HI))
hm_abs = heatmap!(ax_abs0, ζ_xs, ζ_ys, clip.(abs.(W_sheet0));
                  colormap = :viridis, colorrange = (0, U_CAP),
                  nan_color = :transparent)
scatter!(ax_abs0, [0.0], [0.0]; color = :orange, markersize = 10,
         marker = :circle, strokewidth = 1)

ax_phi0 = Axis(fig[1, 2];
               title = "(2a) arg w(ζ), sheet [0]",
               xlabel = "Re ζ", ylabel = "Im ζ",
               aspect = DataAspect(),
               limits = (ζ_RE_LO, ζ_RE_HI, ζ_IM_LO, ζ_IM_HI))
hm_phi = heatmap!(ax_phi0, ζ_xs, ζ_ys, Float64.(angle.(W_sheet0));
                  colormap = :hsv, colorrange = (-π, π),
                  nan_color = :transparent)
scatter!(ax_phi0, [0.0], [0.0]; color = :black, markersize = 10,
         marker = :circle, strokewidth = 1)

ax_z0 = Axis(fig[1, 3];
             title = "(3a) arg u(z), z-plane sheet 0",
             xlabel = "Re z", ylabel = "Im z",
             aspect = DataAspect(),
             limits = (-Z_HALF, Z_HALF, -Z_HALF, Z_HALF))
heatmap!(ax_z0, z_xs, z_ys, Uphi_z_s0;
         colormap = :hsv, colorrange = (-π, π),
         nan_color = :transparent)
lines!(ax_z0, [-Z_HALF, Z_HALF], [0.0, 0.0];
       color = (:black, 0.25), linewidth = 0.3)
lines!(ax_z0, [0.0, 0.0], [-Z_HALF, Z_HALF];
       color = (:black, 0.25), linewidth = 0.3)

# --- Bottom row: alternate sheet --------------------------------------
ax_abs1 = Axis(fig[2, 1];
               title = "(1b) |w(ζ)|, sheet $alt_sheet_label",
               xlabel = "Re ζ", ylabel = "Im ζ",
               aspect = DataAspect(),
               limits = (ζ_RE_LO, ζ_RE_HI, ζ_IM_LO, ζ_IM_HI))
heatmap!(ax_abs1, ζ_xs, ζ_ys, clip.(abs.(W_sheetalt));
         colormap = :viridis, colorrange = (0, U_CAP),
         nan_color = :transparent)
scatter!(ax_abs1, [0.0], [0.0]; color = :orange, markersize = 10,
         marker = :circle, strokewidth = 1)

ax_phi1 = Axis(fig[2, 2];
               title = "(2b) arg w(ζ), sheet $alt_sheet_label",
               xlabel = "Re ζ", ylabel = "Im ζ",
               aspect = DataAspect(),
               limits = (ζ_RE_LO, ζ_RE_HI, ζ_IM_LO, ζ_IM_HI))
heatmap!(ax_phi1, ζ_xs, ζ_ys, Float64.(angle.(W_sheetalt));
         colormap = :hsv, colorrange = (-π, π),
         nan_color = :transparent)
scatter!(ax_phi1, [0.0], [0.0]; color = :black, markersize = 10,
         marker = :circle, strokewidth = 1)

ax_z1 = Axis(fig[2, 3];
             title = "(3b) arg u(z), z-plane sheet $alt_sheet_label",
             xlabel = "Re z", ylabel = "Im z",
             aspect = DataAspect(),
             limits = (-Z_HALF, Z_HALF, -Z_HALF, Z_HALF))
heatmap!(ax_z1, z_xs, z_ys, Uphi_z_alt;
         colormap = :hsv, colorrange = (-π, π),
         nan_color = :transparent)
lines!(ax_z1, [-Z_HALF, Z_HALF], [0.0, 0.0];
       color = (:black, 0.25), linewidth = 0.3)
lines!(ax_z1, [0.0, 0.0], [-Z_HALF, Z_HALF];
       color = (:black, 0.25), linewidth = 0.3)

Colorbar(fig[1, 4], hm_abs; label = "|w|  (capped at $U_CAP)")
Colorbar(fig[2, 4], hm_phi; label = "arg ∈ [-π, π]")

save(OUTPNG, fig)
@printf("\nWrote %s\n", OUTPNG)
