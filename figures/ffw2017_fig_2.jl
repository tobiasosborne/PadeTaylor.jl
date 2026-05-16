# figures/ffw2017_fig_2.jl
#
# Reproduces Fasondini, Fornberg & Weideman 2017, **Figure 2** — the
# headline demonstration of the §2.2 multi-method approach to a
# tronquée P_VI solution.  Three columns, each computing the SAME
# analytic solution by a different method:
#
#   * Column 1 — η-plane (A3 / ADR-less, worklog 041).  PVI in the
#     double-exponential frame `ζ = log z, η = log ζ`.  The
#     branch-point-free region `Re η < log(2π) ≈ 1.838` contains a
#     finite compact rectangle of the solution's domain that maps
#     to "infinitely many ζ-sheets" of the standard transform.  No
#     cut-respecting routing needed inside this region.
#
#   * Column 2 — ζ-plane single sheet, refuse mode (A4 / ADR-0013,
#     worklog 042).  PVI in the standard `z = e^ζ` frame, restricted
#     to one strip `Im ζ ∈ (-π, π]`.  Walker uses
#     `branch_points = (0im,)` with `cross_branch = false`: any wedge
#     candidate whose step crosses the cut from ζ = 0 (along
#     `arg = π`) is rejected.  Demonstrates the "running around the
#     branch point in clockwise and counterclockwise directions"
#     FFW md:178 describes for the bottom-centre frame.
#
#   * Column 3 — ζ-plane multi-sheet, cross mode (A4 + A5, worklog
#     042 + 043).  Same frame as column 2, extended to two strips
#     `Im ζ ∈ (-π, 3π]`.  `cross_branch = true` allows the walker to
#     cross cuts deliberately; `eval_at_sheet` renders the resulting
#     visited tree on a per-pixel sheet basis.  Demonstrates the
#     deliberate crossing FFW md:178 describes for the top-centre
#     frame.
#
# Source: references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
#         FFW2017_painleve_riemann_surfaces_preprint.md:163-189 (the
#         §2.2.2 algorithm); md:195 (Figure 2 caption with parameters,
#         IC, and FFW error estimates).
#
# ## Walker-targets vs render-lattice split (Fig 1 / Fig 6 pattern)
#
# `path_network_solve` walks the walker to EVERY entry in its `grid`
# argument, even if those entries are dense rendering lattice cells.
# For an 80×160 rendering lattice that's ~13k targets — the walker
# spends minutes navigating to each, and the rectangular grid
# transit through pole-dense regions crashes ("all 5 wedge candidates
# failed").
#
# The right pattern (FFW Fig 1 ships this; ffw2017_fig_1.jl §"Manual
# Stage 2"): walk to a SPARSE Stage-1 target set only, then perform
# Stage-2 evaluation manually on the rendering lattice using the
# visited tree's stored per-node Padé approximants.  The package
# accessor `eval_at_sheet` does this for the sheet-aware case;
# `PathNetwork._evaluate_pade` on `sol.visited_pade[k]` does it for
# the sheet-blind case.

using PadeTaylor
using PadeTaylor.SheetTracker:   pVI_transformed_rhs,
                                  pVI_eta_transformed_rhs,
                                  pVI_z_to_η, pVI_η_to_z
using PadeTaylor.CoordTransforms: pV_z_to_ζ, pV_ζ_to_z
using PadeTaylor.PathNetwork:    _evaluate_pade
using CairoMakie
using Printf

# ----------------------------------------------------------------------
# Parameters and ICs (FFW md:195)
# ----------------------------------------------------------------------
const α, β, γ, δ = 4.0, -4.0, 8.0, -8.0
const Z0, U0, UP0 = 10.0+0im, 0.429534600325223+0im, -1.61713114374804e-3+0im

const ORDER       = 30
const ADAPT_TOL   = 1.0e-10
const K_CONS      = 1.0e-3
const MAX_RESC    = 50

# ζ-frame IC (via PV's exponential transform — PVI reuses it; FFW
# md:146).
const ζ0, w0, wp0 = pV_z_to_ζ(Z0, U0, UP0)
# η-frame IC (composition z → ζ → η).
const η0, v0, vp0 = pVI_z_to_η(Z0, U0, UP0)

@printf("FFW 2017 Fig 2: (α, β, γ, δ) = (%.1f, %.1f, %.1f, %.1f)\n",
        α, β, γ, δ)
@printf("                IC (z₀, u, u') = (%s, %s, %s)\n", Z0, U0, UP0)
@printf("                  → (ζ₀, w, w') ≈ (%s, %s, %s)\n", ζ0, w0, wp0)
@printf("                  → (η₀, v, v') ≈ (%s, %s, %s)\n", η0, v0, vp0)
@printf("                Re η₀ = %.3f  <  log(2π) = %.3f  ✓\n",
        real(η0), log(2π))

# ----------------------------------------------------------------------
# Window definitions
# ----------------------------------------------------------------------

# η-plane window — inside the branch-point-free region.
# v1 corner: FFW Fig 2 shows TWO strips (Im η ∈ (-π, 3π]), but our
# uniform-target walker struggles with the high-|Im η| × small-Re η
# corner under the dense PVI pole structure.  We render a SINGLE
# strip Im η ∈ (-π, π] (the principal sheet in the folded η
# representation) for v1; widening the window is the same kind of
# Poisson-disk-node lift as bead `padetaylor-zwh` (FFW Fig 1).
const η_RE_LO, η_RE_HI = -0.5, log(2π) - 0.05
const η_IM_LO, η_IM_HI = -π + 0.05, π - 0.05

# ζ-plane window for Columns 2 + 3.  Re ζ runs from 0 (avoid crossing
# the negative-real cut by construction: the rendering lattice stays
# on the positive-real side of the origin branch, so refuse mode is
# active but never triggers).  Im ζ window: one strip for col 2, two
# strips for col 3.
const ζ_RE_LO, ζ_RE_HI = 0.05, 3.0
const ζ2_IM_LO, ζ2_IM_HI = -π + 0.05, π - 0.05
const ζ3_IM_LO, ζ3_IM_HI = -π + 0.05, 3π - 0.05

# Stage-2 evaluation lattice density.
const NX, NY_SINGLE = 60, 60
const NY_DOUBLE     = 120
const NY_η          = NY_SINGLE

# ----------------------------------------------------------------------
# Build SPARSE rectangular Stage-1 target grids (NOT the render
# lattice).  ~150-300 targets each — the walker navigates these
# robustly under adaptive_ffw + refuse / cross mode.
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

# ----------------------------------------------------------------------
# Solve A — η-plane (column 1)
# ----------------------------------------------------------------------
println("\n== Solve A: η-plane (A3) ==")
const η_targets = rect_targets(η_RE_LO, η_RE_HI, η_IM_LO, η_IM_HI, 0.35, 0.7)
@printf("  Stage-1 targets: %d\n", length(η_targets))

const f_η = pVI_eta_transformed_rhs(α, β, γ, δ)
const prob_η = PadeTaylorProblem(f_η, (v0, vp0),
                                 (η0, complex(η_RE_HI, η_IM_HI));
                                 order = ORDER)
t0_η = time()
const sol_η = path_network_solve(prob_η, η_targets;
                                  h = 0.3,
                                  step_size_policy = :adaptive_ffw,
                                  adaptive_tol = ADAPT_TOL,
                                  k_conservative = K_CONS,
                                  max_rescales = MAX_RESC,
                                  max_steps_per_target = 3000)
@printf("  Stage 1 in %.2f s; %d visited nodes\n",
        time() - t0_η, length(sol_η.visited_z))

# ----------------------------------------------------------------------
# Solve B — ζ-plane single sheet, refuse mode (column 2)
# ----------------------------------------------------------------------
println("\n== Solve B: ζ-plane refuse mode (A4) ==")
const ζ2_targets = rect_targets(ζ_RE_LO, ζ_RE_HI, ζ2_IM_LO, ζ2_IM_HI, 0.35, 0.7)
@printf("  Stage-1 targets: %d\n", length(ζ2_targets))

const f_ζ = pVI_transformed_rhs(α, β, γ, δ)
const prob_ζ2 = PadeTaylorProblem(f_ζ, (w0, wp0),
                                  (ζ0, complex(ζ_RE_HI, ζ2_IM_HI));
                                  order = ORDER)
t0_ζ2 = time()
const sol_ζ2 = path_network_solve(prob_ζ2, ζ2_targets;
                                   h = 0.3,
                                   step_size_policy = :adaptive_ffw,
                                   adaptive_tol = ADAPT_TOL,
                                   k_conservative = K_CONS,
                                   max_rescales = MAX_RESC,
                                   max_steps_per_target = 3000,
                                   branch_points     = (0.0 + 0.0im,),
                                   branch_cut_angles = π)
@printf("  Stage 1 in %.2f s; %d visited nodes\n",
        time() - t0_ζ2, length(sol_ζ2.visited_z))

# ----------------------------------------------------------------------
# Solve C — ζ-plane multi-sheet, cross mode (column 3)
# ----------------------------------------------------------------------
println("\n== Solve C: ζ-plane cross mode (A4 + A5) ==")
# Rectangular targets in the positive-Re ζ region (avoids cornering)
# PLUS a winding ring around ζ = 0 in the NEGATIVE-Re region that
# straddles the cut.  A target in upper-half then lower-half (or
# vice versa) at the same negative Re forces the walker to step
# from positive Im to negative Im at Re < 0 — crossing the cut and
# bumping the sheet counter.
const ζ3_targets_base = rect_targets(ζ_RE_LO, ζ_RE_HI, ζ3_IM_LO, ζ3_IM_HI,
                                      0.35, 0.7)
const winding_ring = ComplexF64[
    -0.5 + 0.5im, -0.5 - 0.5im,
    -1.0 + 0.5im, -1.0 - 0.5im,
    -1.0 + 1.0im, -1.0 - 1.0im,
    -0.5 + 1.0im, -0.5 - 1.0im,
]
const ζ3_targets = vcat(ζ3_targets_base, winding_ring)
@printf("  Stage-1 targets: %d (incl. %d winding ring)\n",
        length(ζ3_targets), length(winding_ring))

const prob_ζ3 = PadeTaylorProblem(f_ζ, (w0, wp0),
                                  (ζ0, complex(ζ_RE_HI, ζ3_IM_HI));
                                  order = ORDER)
t0_ζ3 = time()
const sol_ζ3 = path_network_solve(prob_ζ3, ζ3_targets;
                                   h = 0.3,
                                   step_size_policy = :adaptive_ffw,
                                   adaptive_tol = ADAPT_TOL,
                                   k_conservative = K_CONS,
                                   max_rescales = MAX_RESC,
                                   max_steps_per_target = 3000,
                                   branch_points     = (0.0 + 0.0im,),
                                   branch_cut_angles = π,
                                   cross_branch      = true)
@printf("  Stage 1 in %.2f s; %d visited nodes\n",
        time() - t0_ζ3, length(sol_ζ3.visited_z))

const sheet_counts = Dict{Vector{Int},Int}()
for s in sol_ζ3.visited_sheet
    sheet_counts[s] = get(sheet_counts, s, 0) + 1
end
println("  Visited-sheet distribution: ", sheet_counts)

# ----------------------------------------------------------------------
# Manual Stage 2 for sheet-blind solves: nearest-visited Padé eval.
# ----------------------------------------------------------------------
function stage2_eval_blind(sol::PathNetworkSolution, points)
    nv = length(sol.visited_z)
    out = Vector{ComplexF64}(undef, length(points))
    for (i, zf) in enumerate(points)
        best_k = 1
        best_d = abs(zf - sol.visited_z[1])
        for k in 2:nv
            d = abs(zf - sol.visited_z[k])
            if d < best_d
                best_d = d
                best_k = k
            end
        end
        hv = sol.visited_h[best_k]
        if best_d > hv
            out[i] = complex(NaN, NaN)
        else
            t = (zf - sol.visited_z[best_k]) / hv
            out[i] = _evaluate_pade(sol.visited_pade[best_k], t)
        end
    end
    return out
end

# ----------------------------------------------------------------------
# Render lattices
# ----------------------------------------------------------------------
const η_xs  = collect(range(η_RE_LO, η_RE_HI; length = NX))
const η_ys  = collect(range(η_IM_LO, η_IM_HI; length = NY_η))
const η_lattice = ComplexF64[complex(x, y) for y in η_ys for x in η_xs]

const ζ_xs  = collect(range(ζ_RE_LO, ζ_RE_HI; length = NX))
const ζ2_ys = collect(range(ζ2_IM_LO, ζ2_IM_HI; length = NY_SINGLE))
const ζ2_lattice = ComplexF64[complex(x, y) for y in ζ2_ys for x in ζ_xs]

const ζ3_ys = collect(range(ζ3_IM_LO, ζ3_IM_HI; length = NY_DOUBLE))

t0_eval = time()
const W_η  = reshape(stage2_eval_blind(sol_η,  η_lattice),  NX, NY_η)
const W_ζ2 = reshape(stage2_eval_blind(sol_ζ2, ζ2_lattice), NX, NY_SINGLE)

# Column 3: per-sheet eval via eval_at_sheet on sol_ζ3.  The
# "alternate sheet" populated by the walker depends on which way the
# winding ring's natural traversal order winds — could be sheet [+1]
# (CCW) or [-1] (CW).  We pick whichever has the larger node count.
const alt_sheet_label = let
    counts = Dict{Vector{Int},Int}()
    for s in sol_ζ3.visited_sheet
        counts[s] = get(counts, s, 0) + 1
    end
    # Among non-[0] sheets, the one with the most visited nodes.
    nonzero = [(s, n) for (s, n) in counts if s != [0]]
    if isempty(nonzero)
        Int[]
    else
        sort!(nonzero; by = x -> -x[2])
        first(nonzero)[1]
    end
end
@printf("  Alternate-sheet rendering target: %s\n", alt_sheet_label)

W_ζ3_sheet0   = Matrix{ComplexF64}(undef, NX, NY_DOUBLE)
W_ζ3_sheetalt = Matrix{ComplexF64}(undef, NX, NY_DOUBLE)
for j in 1:NY_DOUBLE, i in 1:NX
    z = complex(ζ_xs[i], ζ3_ys[j])
    u0_at, _ = eval_at_sheet(sol_ζ3, z, [0])
    W_ζ3_sheet0[i, j]   = u0_at
    if !isempty(alt_sheet_label)
        u_alt, _ = eval_at_sheet(sol_ζ3, z, alt_sheet_label)
        W_ζ3_sheetalt[i, j] = u_alt
    else
        W_ζ3_sheetalt[i, j] = complex(NaN, NaN)
    end
end
@printf("  Stage 2 eval in %.2f s\n", time() - t0_eval)

const cov_η  = count(isfinite, abs.(W_η))
const cov_ζ2 = count(isfinite, abs.(W_ζ2))
const cov_s0  = count(isfinite, abs.(W_ζ3_sheet0))
const cov_alt = count(isfinite, abs.(W_ζ3_sheetalt))
@printf("  Coverage:  η: %d/%d (%.1f%%), ζ-refuse: %d/%d (%.1f%%), ζ-sheet0: %d/%d (%.1f%%), ζ-sheet%s: %d/%d (%.1f%%)\n",
        cov_η, length(W_η), 100*cov_η/length(W_η),
        cov_ζ2, length(W_ζ2), 100*cov_ζ2/length(W_ζ2),
        cov_s0, length(W_ζ3_sheet0), 100*cov_s0/length(W_ζ3_sheet0),
        alt_sheet_label, cov_alt, length(W_ζ3_sheetalt),
        100*cov_alt/length(W_ζ3_sheetalt))

# ----------------------------------------------------------------------
# Render — three columns, |u| heatmap in each.
# ----------------------------------------------------------------------
const U_CAP = 6.0
clip(x) = isfinite(x) ? min(x, U_CAP) : NaN

const OUTPNG = joinpath(@__DIR__, "output", "ffw2017_fig_2.png")
mkpath(dirname(OUTPNG))

fig = Figure(size = (1800, 1100))

ax_η = Axis(fig[1, 1];
            title = "(1) η-plane (A3): |v(η)|",
            xlabel = "Re η", ylabel = "Im η",
            aspect = DataAspect())
hm_η = heatmap!(ax_η, η_xs, η_ys, clip.(abs.(W_η));
                colormap = :viridis, colorrange = (0, U_CAP),
                nan_color = :transparent)
vlines!(ax_η, [log(2π)]; color = :red, linestyle = :dash)
scatter!(ax_η, [real(η0)], [imag(η0)]; color = :white,
         markersize = 12, marker = :star5, strokewidth = 1)
text!(ax_η, real(η0) + 0.05, imag(η0); text = "η₀", color = :white,
      fontsize = 11)

ax_ζ2 = Axis(fig[1, 2];
             title = "(2) ζ-plane sheet 0, refuse (A4): |w(ζ)|",
             xlabel = "Re ζ", ylabel = "Im ζ",
             aspect = DataAspect())
hm_ζ2 = heatmap!(ax_ζ2, ζ_xs, ζ2_ys, clip.(abs.(W_ζ2));
                  colormap = :viridis, colorrange = (0, U_CAP),
                  nan_color = :transparent)
scatter!(ax_ζ2, [0.0], [0.0]; color = :orange, markersize = 12,
         marker = :circle, strokewidth = 1)
scatter!(ax_ζ2, [real(ζ0)], [imag(ζ0)]; color = :white,
         markersize = 12, marker = :star5, strokewidth = 1)
text!(ax_ζ2, real(ζ0) + 0.05, imag(ζ0); text = "ζ₀", color = :white,
      fontsize = 11)

ax_ζ3a = Axis(fig[1, 3];
              title = "(3a) ζ sheet [0], cross (A4+A5)",
              xlabel = "Re ζ", ylabel = "Im ζ",
              aspect = DataAspect())
heatmap!(ax_ζ3a, ζ_xs, ζ3_ys, clip.(abs.(W_ζ3_sheet0));
          colormap = :viridis, colorrange = (0, U_CAP),
          nan_color = :transparent)
scatter!(ax_ζ3a, [0.0], [0.0]; color = :orange, markersize = 12,
         marker = :circle, strokewidth = 1)

ax_ζ3b = Axis(fig[2, 3];
              title = "(3b) ζ sheet $alt_sheet_label, cross (A4+A5)",
              xlabel = "Re ζ", ylabel = "Im ζ",
              aspect = DataAspect())
heatmap!(ax_ζ3b, ζ_xs, ζ3_ys, clip.(abs.(W_ζ3_sheetalt));
          colormap = :viridis, colorrange = (0, U_CAP),
          nan_color = :transparent)
scatter!(ax_ζ3b, [0.0], [0.0]; color = :orange, markersize = 12,
         marker = :circle, strokewidth = 1)

Label(fig[2, 1:2], """
FFW 2017 Fig 2 — three-method reproduction of a tronquée P_VI solution.
(α,β,γ,δ) = $((α, β, γ, δ));  z-plane IC u(10) = $U0, u'(10) = $UP0.
Method A (η-plane, A3):       $(length(sol_η.visited_z)) visited nodes, $(round(time()-t0_eval; digits=1))s eval.
Method B (ζ-plane refuse, A4): $(length(sol_ζ2.visited_z)) visited nodes.
Method C (ζ-plane cross + A5): $(length(sol_ζ3.visited_z)) visited nodes;
   visited_sheet distribution: $(sheet_counts).
   Sheet-[0] coverage $(round(100*cov_s0/length(W_ζ3_sheet0); digits=1))%;
   sheet-$alt_sheet_label coverage $(round(100*cov_alt/length(W_ζ3_sheetalt); digits=1))%.""";
      tellwidth = false)

Colorbar(fig[1, 4], hm_η; label = "|u|  (capped at $U_CAP)")

save(OUTPNG, fig)
@printf("\nWrote %s\n", OUTPNG)
