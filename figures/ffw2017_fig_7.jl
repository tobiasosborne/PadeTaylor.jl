# figures/ffw2017_fig_7.jl
#
# Reproduces Fasondini, Fornberg & Weideman 2017, **Figure 7** — a
# generic P_VI solution in three coordinate frames simultaneously:
# η, ζ and z.  This is the last (T5) figure in the FFW arc.
#
# Source: references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
#         FFW2017_painleve_riemann_surfaces_preprint.md:281-301.  The
#         caption verbatim (md:301):
#
# > Figure 7: A generic PVI solution in the η, ζ and z planes
# > (ζ = e^η, z = e^ζ).  The solution has parameters
# > (α, β, γ, δ) = (1, −1, 3/4, −3/2) and ICs z₀ = 2, u(z₀) = 3/2 and
# > u′(2) = −1.  The error estimates for the solution on the strips
# > indicated in the η-plane are, from bottom to top: 5e-8, 4e-4 and
# > 8e-4.
#
# The supporting paragraph at md:281 records WHY the z-plane uses
# `log10|u|` instead of the `|u|` employed in η and ζ:
#
# > It follows from Table 3 that for a P_VI solution with a pole at
# > z_0, where |z_0| > 1, the residue of the pole in the z-plane is
# > larger than the corresponding residues in the η and ζ planes by a
# > factor of at least |z_0|.  We therefore found it necessary to plot
# > log_10 |u| in the z-plane in Figure 7 (column 3), instead of |u|,
# > which is what we plot in the η and ζ planes (columns 1 and 2,
# > respectively).
#
# ## What Fig 7 reveals (vs Figs 2 + 3)
#
# Figures 2 and 3 show a *tronquée* PVI solution (large pole-free
# sector).  Figure 7 shows the GENERIC case — no large pole-free
# sector; instead dense oblique pole fields throughout the rendered
# windows.  This means:
#
#   * In the ζ-plane the pole density is high everywhere, NOT just
#     close to the branch points; the `node_separation = R(ζ)`
#     non-uniform Stage-1 placement (ADR-0012) is REQUIRED for the
#     walker to handle the dense pole field without bogging.  (Figs 2
#     and 3 dodged this with their sparse rectangular base + winding
#     ring because the tronquée's pole-free sector covers most of the
#     window.)
#
#   * In the z-plane the residue magnification `|z_0|` (Table 3) at
#     poles `z_0` with `|z_0| > 1` means the rendered `|u(z)|` has
#     dynamic range spanning many orders of magnitude even after the
#     usual viridis cap.  FFW switch to `log_10|u|` for column 3 so
#     pole spikes and the bulk solution coexist on one heatmap.
#
# ## Two solves
#
#   * Solve E (column 1) — η-plane via `pVI_eta_transformed_rhs`.  We
#     walk inside the branch-point-free region `Re η < log(2π)`
#     (FFW md:157).  IC `(η₀, v₀, vp₀) ≈ (-0.367, 1.5, -1.386)` is
#     well inside this region.
#
#   * Solve F (columns 2 + 3) — ζ-plane via `pVI_transformed_rhs` with
#     `cross_branch = true` and a winding ring.  Column 2 reads the
#     resulting `W` array directly; column 3 is derived by Cartesian
#     resample of `W` onto a z-plane lattice via `z = e^ζ` and
#     `log10|u|` colour mapping.
#
# ## Branch-point inventory (NEW v1 risk handled inline)
#
# PVI has TWO fixed branch points: `z = 0` (mapped to `ζ = -∞`, away)
# AND `z = 1` (mapped to `ζ = 2πi·k`, k ∈ ℤ — a lattice that survives
# in the ζ-plane; FFW md:141).  For the tronquée Figs 2 + 3 we shipped
# `branch_points = (0.0+0.0im,)` only, because the pole-free sector
# kept the walker on the LOW-|Im ζ| portion of the window where the
# `ζ = 2πi` lattice point is irrelevant.  For the GENERIC case with
# `Im ζ ∈ (-π, 3π]` the walker DOES approach `ζ = 2πi`, so we add
# that lattice point (and its mirror `ζ = -2πi`) explicitly to
# `branch_points`.  All three cuts run along `arg = π` (Julia's `log`
# default cut convention rotated to each branch point).  See
# `src/BranchTracker.jl` lines 71-81 — the kwarg expects a tuple of
# `Complex` of arbitrary length with `branch_cut_angles` matching.
#
# The walker's `visited_sheet` now tracks a length-3 `Vector{Int}`
# (one counter per branch), but for the purposes of the
# `eval_at_sheet` query in the alternate-sheet panel we look at the
# FIRST counter only (the `ζ = 0` branch — same as Fig 3) because
# the load-bearing claim is the original `ζ = 0` cut-crossing.  Sheet
# inspections are reported in the figure-script `@printf` output and
# in the worklog.
#
# ## rng_seed sensitivity (honest v1 corner)
#
# The walker's `shuffle(rng, targets)` step (FW 2011 line 156) means
# walker outcomes depend on `rng_seed`.  For Fig 3 a small sweep of
# `{0, 1, 2, 17, 42}` produced different alternate-sheet populations.
# We start Fig 7 at `rng_seed = 0` (the default) and verify the walker
# successfully produces sheet-[1] visited nodes; if it does not, we
# fall back to `rng_seed = 2` (the Fig 3 choice).  The chosen seed
# and the sheet distribution are reported in the script output and
# pinned in the test.

using PadeTaylor
using PadeTaylor.SheetTracker:   pVI_transformed_rhs,
                                  pVI_eta_transformed_rhs,
                                  pVI_z_to_η, pVI_η_to_z
using PadeTaylor.CoordTransforms: pV_z_to_ζ, pV_ζ_to_z
using CairoMakie
using Printf

# ----------------------------------------------------------------------
# Parameters and ICs (FFW md:301 verbatim)
# ----------------------------------------------------------------------
const α, β, γ, δ = 1.0, -1.0, 0.75, -1.5
const Z0, U0, UP0 = 2.0 + 0.0im, 1.5 + 0.0im, -1.0 + 0.0im

const ORDER       = 30
const ADAPT_TOL   = 1.0e-10
const K_CONS      = 1.0e-3
const MAX_RESC    = 50

# ζ-frame IC via the PV exp transform (PVI inherits this; FFW md:146).
const ζ0, w0, wp0 = pV_z_to_ζ(Z0, U0, UP0)
# η-frame IC via the double-exp composition (FFW md:151).
const η0, v0, vp0 = pVI_z_to_η(Z0, U0, UP0)

@printf("FFW 2017 Fig 7: (α, β, γ, δ) = (%.2f, %.2f, %.2f, %.2f)\n",
        α, β, γ, δ)
@printf("                IC (z₀, u, u') = (%s, %s, %s)\n", Z0, U0, UP0)
@printf("                  → (ζ₀, w, w') ≈ (%s, %s, %s)\n", ζ0, w0, wp0)
@printf("                  → (η₀, v, v') ≈ (%s, %s, %s)\n", η0, v0, vp0)
@printf("                Re η₀ = %.3f < log(2π) = %.3f  ✓\n",
        real(η0), log(2π))

# ----------------------------------------------------------------------
# Windows
# ----------------------------------------------------------------------

# η-window — strictly inside `Re η < log(2π)`.
const η_RE_LO, η_RE_HI = -1.0, log(2π) - 0.05
const η_IM_LO, η_IM_HI = -π + 0.05, π - 0.05

# ζ-window for columns 2 + 3.  Two strips Im ζ ∈ (-π, 3π] (sheet [0]
# below the cut and sheet [1] above) — same shape as Fig 3.  Re ζ
# bounded a bit tighter than Fig 3 (`[-0.8, 1.5]`) because (a) the
# generic case has dense poles everywhere so the walker doesn't need
# extra Re ζ to reveal structure and (b) the z-plane half-window is
# `exp(Re ζ_HI) ≈ 4.5` which already contains the IC at `z = 2`.
const ζ_RE_LO, ζ_RE_HI = -0.8, 1.5
const ζ_IM_LO, ζ_IM_HI = -π + 0.05, 3π - 0.05

# Stage-2 render-lattice density.  ζ-plane gets a tighter mesh than
# Figs 2/3 because the dense poles need finer rendering.
const NX_η, NY_η = 90, 110
const NX_ζ, NY_ζ = 90, 220

# z-plane render half-window: contains `|z| ≤ exp(Re ζ_HI) ≈ 4.5`.
const Z_RMAX = exp(ζ_RE_HI)
const Z_HALF = 1.05 * Z_RMAX
const Z_NX, Z_NY = 240, 240

# Non-uniform Stage-1 spacing for the ζ-frame walker (ADR-0012; Fig 6
# style).  Same shape `R(ζ) = max(0.05, (4 - Re ζ)/10)` as Fig 6 — at
# `Re ζ = -0.8` gives `R = 0.48`, at `Re ζ = 1.5` gives `R = 0.25`.
# Generic PVI doesn't have Fig 6's smooth-region growth so a constant
# floor + linear slope is sufficient.
R_ζ(ζ) = max(0.05, (4.0 - real(ζ)) / 10.0)

# ----------------------------------------------------------------------
# Solve E — η-plane (column 1)
# ----------------------------------------------------------------------
println("\n== Solve E: η-plane (A3) ==")

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

const η_targets = rect_targets(η_RE_LO, η_RE_HI, η_IM_LO, η_IM_HI,
                                0.35, 0.6)
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
                                  max_steps_per_target = 4000)
@printf("  Stage 1 in %.2f s; %d visited nodes\n",
        time() - t0_η, length(sol_η.visited_z))

# ----------------------------------------------------------------------
# Solve F — ζ-plane multi-sheet (columns 2 + 3)
# ----------------------------------------------------------------------
println("\n== Solve F: ζ-plane cross mode + non-uniform nodes (A1+A2+A4+A5) ==")

# Stage-1 targets: R(ζ)-spaced rectangular base (positive-Re region
# only — the cut at Re < 0 is crossed via the winding ring) + 8-target
# winding ring straddling the cut.  Same architecture as Fig 3 but
# with R(ζ)-spaced base (Fig 6 pattern) for the dense generic case.
function build_R_targets(re_lo, re_hi, im_lo, im_hi)
    tgs = ComplexF64[]
    re = re_lo
    while re ≤ re_hi + 1e-12
        rh = R_ζ(re + 0.0im)
        imv = im_lo
        while imv ≤ im_hi + 1e-12
            push!(tgs, complex(re, imv))
            imv += rh
        end
        re += rh
    end
    return tgs
end

const targets_base = build_R_targets(0.05, ζ_RE_HI, ζ_IM_LO, ζ_IM_HI)
# Winding ring identical to Fig 3 — proven cut-crossing geometry.
const winding_ring = ComplexF64[
    -0.5 + 0.5im, -0.5 - 0.5im,
    -1.0 + 0.5im, -1.0 - 0.5im,
    -1.0 + 1.0im, -1.0 - 1.0im,
    -0.5 + 1.0im, -0.5 - 1.0im,
]
const stage1_targets = vcat(targets_base, winding_ring)
@printf("  Stage-1 targets: %d (incl. %d winding ring, R(ζ)-spaced)\n",
        length(stage1_targets), length(winding_ring))

# Branch-point inventory.  Try three lattice points (k = -1, 0, +1).
# All cuts run along arg = π (Julia log convention).  See module
# docstring; src/BranchTracker.jl confirms tuples of arbitrary length
# are supported (the `length(branches) == length(cut_angles)` assertion
# in `any_cut_crossed` enforces matching length).
const BRANCHES = (0.0 + 0.0im, 0.0 + 2π*im, 0.0 - 2π*im)
const CUT_ANGLES = (Float64(π), Float64(π), Float64(π))

const f_ζ = pVI_transformed_rhs(α, β, γ, δ)
const prob_ζ = PadeTaylorProblem(f_ζ, (w0, wp0),
                                  (ζ0, complex(ζ_RE_HI, ζ_IM_HI));
                                  order = ORDER)

t0_ζ = time()
const RNG_SEED = 1
const sol_ζ = path_network_solve(prob_ζ, stage1_targets;
                                  h = R_ζ(ζ0),
                                  node_separation = R_ζ,
                                  step_size_policy = :adaptive_ffw,
                                  adaptive_tol = ADAPT_TOL,
                                  k_conservative = K_CONS,
                                  max_rescales = MAX_RESC,
                                  max_steps_per_target = 4000,
                                  branch_points     = BRANCHES,
                                  branch_cut_angles = CUT_ANGLES,
                                  cross_branch      = true,
                                  rng_seed          = RNG_SEED)
@printf("  Stage 1 in %.2f s; %d visited nodes\n",
        time() - t0_ζ, length(sol_ζ.visited_z))

const sheet_counts = Dict{Vector{Int},Int}()
for s in sol_ζ.visited_sheet
    sheet_counts[s] = get(sheet_counts, s, 0) + 1
end
println("  Visited-sheet distribution (per-branch [ζ=0, ζ=2πi, ζ=-2πi]): ",
        sheet_counts)

# ----------------------------------------------------------------------
# Pick the alternate sheet: the non-`zeros` sheet vector with the most
# visited nodes.  Same heuristic as Figs 2 + 3 but generalised for the
# length-3 sheet vector.
# ----------------------------------------------------------------------
const SHEET_ZERO = zeros(Int, length(BRANCHES))
const alt_sheet_label = let
    nonzero = [(s, n) for (s, n) in sheet_counts if s != SHEET_ZERO]
    isempty(nonzero) ? Int[] :
        first(sort!(nonzero; by = x -> -x[2]))[1]
end
@printf("  Alternate-sheet rendering target: %s\n", alt_sheet_label)

# ----------------------------------------------------------------------
# Stage-2 evaluation lattices
# ----------------------------------------------------------------------
const η_xs = collect(range(η_RE_LO, η_RE_HI; length = NX_η))
const η_ys = collect(range(η_IM_LO, η_IM_HI; length = NY_η))

const ζ_xs = collect(range(ζ_RE_LO, ζ_RE_HI; length = NX_ζ))
const ζ_ys = collect(range(ζ_IM_LO, ζ_IM_HI; length = NY_ζ))

# η-plane Stage-2 eval (sheet-blind via the sol_η accessor).
println("\n== Stage 2: η + ζ Stage-2 eval ==")
V_η = Matrix{ComplexF64}(undef, NX_η, NY_η)
t0_e = time()
for j in 1:NY_η, i in 1:NX_η
    z = complex(η_xs[i], η_ys[j])
    u, _ = eval_at(sol_η, z; extrapolate = true)
    V_η[i, j] = u
end
@printf("  η Stage-2 eval in %.2f s\n", time() - t0_e)

# ζ-plane Stage-2 eval on TWO sheets.
W_sheet0 = Matrix{ComplexF64}(undef, NX_ζ, NY_ζ)
W_sheetalt = Matrix{ComplexF64}(undef, NX_ζ, NY_ζ)
t0_e = time()
for j in 1:NY_ζ, i in 1:NX_ζ
    z = complex(ζ_xs[i], ζ_ys[j])
    u0_at, _ = eval_at_sheet(sol_ζ, z, SHEET_ZERO; extrapolate = true)
    W_sheet0[i, j] = u0_at
    if !isempty(alt_sheet_label)
        ualt, _ = eval_at_sheet(sol_ζ, z, alt_sheet_label;
                                extrapolate = true)
        W_sheetalt[i, j] = ualt
    else
        W_sheetalt[i, j] = complex(NaN, NaN)
    end
end
@printf("  ζ Stage-2 eval in %.2f s\n", time() - t0_e)

const cov_η   = count(isfinite, abs.(V_η))
const cov_s0  = count(isfinite, abs.(W_sheet0))
const cov_alt = count(isfinite, abs.(W_sheetalt))
@printf("  Coverage:  η %d/%d (%.1f%%), ζ sheet0 %d/%d (%.1f%%), ζ alt %d/%d (%.1f%%)\n",
        cov_η, length(V_η), 100 * cov_η / length(V_η),
        cov_s0, length(W_sheet0), 100 * cov_s0 / length(W_sheet0),
        cov_alt, length(W_sheetalt), 100 * cov_alt / length(W_sheetalt))

# ----------------------------------------------------------------------
# Render helpers — `bilinear_w` (bilinear interp of `w(ζ)` from the
# ζ-grid) and `zplane_log_modulus` (the load-bearing z-plane rendering
# invariant, FFW md:281).  Extracted into a Makie-free kernel so the
# test file can `include` them too and exercise the SAME code path
# (FF7.1.8 mutation-prove M5).  See `_ffw2017_fig_7_helpers.jl`.
# ----------------------------------------------------------------------
include(joinpath(@__DIR__, "_ffw2017_fig_7_helpers.jl"))

# Aggregate the z-plane log10|u| heatmap for one sheet.  Uses the
# helper so the test's per-cell pin and the figure render share code.
const z_xs = collect(range(-Z_HALF, Z_HALF; length = Z_NX))
const z_ys = collect(range(-Z_HALF, Z_HALF; length = Z_NY))

function render_zplane_log(W::Matrix{ComplexF64}, sheet_lift::Integer)
    out = fill(NaN, Z_NX, Z_NY)
    for jx in 1:Z_NY, ix in 1:Z_NX
        out[ix, jx] = zplane_log_modulus(W, complex(z_xs[ix], z_ys[jx]),
                                          sheet_lift,
                                          ζ_xs, ζ_ys,
                                          ζ_RE_LO, ζ_RE_HI)
    end
    return out
end

# Choose `sheet_lift` per panel.  Sheet 0 → lift 0.  Alternate sheet:
# inspect the FIRST branch-counter (the ζ = 0 branch) and lift the
# z-plane arg accordingly (mirroring Fig 3's choice).
const alt_lift = isempty(alt_sheet_label) ? 0 :
                  (first(alt_sheet_label) == 0 ? 1 : Int(first(alt_sheet_label)))

t0_z = time()
const Ulog_z_s0  = render_zplane_log(W_sheet0,   0)
const Ulog_z_alt = render_zplane_log(W_sheetalt, alt_lift)
@printf("  z-plane log10|u| render in %.2f s\n", time() - t0_z)

# ----------------------------------------------------------------------
# Per-strip symmetry-based error estimates (advisory, NOT load-bearing
# per the bead's "do NOT try to beat the 5e-8 / 4e-4 / 8e-4 pins").
# FFW md:301 uses the Schwarz-reflection comparison on the η-frame
# bottom-to-top strips.  We restrict to the rendered single strip
# `Im η ∈ (-π, π]` so this is one number, not three.
# ----------------------------------------------------------------------
let
    err_max  = 0.0
    err_med  = 0.0
    samples  = ComplexF64[]
    for re in range(η_RE_LO + 0.1, η_RE_HI - 0.05; length = 10)
        for im in range(0.1, η_IM_HI - 0.1; length = 4)
            push!(samples, complex(re, im))
        end
    end
    diffs = Float64[]
    for s in samples
        u_top,    _ = eval_at(sol_η, s; extrapolate = true)
        u_bottom, _ = eval_at(sol_η, conj(s); extrapolate = true)
        isfinite(real(u_top)) && isfinite(real(u_bottom)) || continue
        push!(diffs, abs(u_top - conj(u_bottom)))
    end
    if !isempty(diffs)
        err_max = maximum(diffs)
        err_med = sort(diffs)[div(length(diffs)+1, 2)]
    end
    @info "η-plane Schwarz-reflection error (max / median over $(length(diffs)) samples)" err_max err_med
end

# ----------------------------------------------------------------------
# Render — three columns × two rows (sheet 0 top, alt sheet bottom).
#
#   Col 1: |v(η)| modulus (sheet 0 only, alt-row blank)
#   Col 2: |w(ζ)| modulus per sheet
#   Col 3: log10|u(z)| per sheet
#
# The η-plane is a single rectangular face (FFW Fig 7 col 1 is one
# composite cell containing three bottom-to-top strips with error
# pins — we ship the principal-sheet face only as v1 corner; see
# worklog 047).
# ----------------------------------------------------------------------
const ABS_CAP = 6.0
const LOG_LO, LOG_HI = -1.5, 2.0
clip(x) = isfinite(x) ? min(x, ABS_CAP) : NaN

const OUTPNG = joinpath(@__DIR__, "output", "ffw2017_fig_7.png")
mkpath(dirname(OUTPNG))

fig = Figure(size = (1700, 1050))
Label(fig[0, 1:3],
      "FFW 2017 Fig 7 — generic PVI on η, ζ and z planes, " *
      "(α,β,γ,δ) = $((α, β, γ, δ)),  IC (z,u,u') = $((Z0, U0, UP0))";
      fontsize = 15, padding = (0, 0, 4, 4))

# ----- Top row: principal sheet ---------------------------------------
ax_η = Axis(fig[1, 1];
            title = "(1) |v(η)|, η-plane (Re η < log 2π)",
            xlabel = "Re η", ylabel = "Im η",
            aspect = DataAspect(),
            limits = (η_RE_LO, η_RE_HI, η_IM_LO, η_IM_HI))
hm_η = heatmap!(ax_η, η_xs, η_ys, clip.(abs.(V_η));
                colormap = :viridis, colorrange = (0, ABS_CAP),
                nan_color = :transparent)
vlines!(ax_η, [log(2π)]; color = :red, linestyle = :dash)
scatter!(ax_η, [real(η0)], [imag(η0)]; color = :white,
         markersize = 11, marker = :star5, strokewidth = 1)

ax_ζ0 = Axis(fig[1, 2];
             title = "(2a) |w(ζ)|, ζ-plane sheet 0",
             xlabel = "Re ζ", ylabel = "Im ζ",
             aspect = DataAspect(),
             limits = (ζ_RE_LO, ζ_RE_HI, ζ_IM_LO, ζ_IM_HI))
hm_ζ = heatmap!(ax_ζ0, ζ_xs, ζ_ys, clip.(abs.(W_sheet0));
                colormap = :viridis, colorrange = (0, ABS_CAP),
                nan_color = :transparent)
scatter!(ax_ζ0, [0.0], [0.0]; color = :orange, markersize = 10,
         marker = :circle, strokewidth = 1)
scatter!(ax_ζ0, [0.0], [2π]; color = :orange, markersize = 10,
         marker = :circle, strokewidth = 1)

ax_z0 = Axis(fig[1, 3];
             title = "(3a) log₁₀|u(z)|, z-plane sheet 0",
             xlabel = "Re z", ylabel = "Im z",
             aspect = DataAspect(),
             limits = (-Z_HALF, Z_HALF, -Z_HALF, Z_HALF))
hm_z = heatmap!(ax_z0, z_xs, z_ys, Ulog_z_s0;
                colormap = :viridis, colorrange = (LOG_LO, LOG_HI),
                nan_color = :transparent)
lines!(ax_z0, [-Z_HALF, Z_HALF], [0.0, 0.0];
       color = (:black, 0.25), linewidth = 0.3)
lines!(ax_z0, [0.0, 0.0], [-Z_HALF, Z_HALF];
       color = (:black, 0.25), linewidth = 0.3)
scatter!(ax_z0, [real(Z0)], [imag(Z0)]; color = :white,
         markersize = 11, marker = :star5, strokewidth = 1)

# ----- Bottom row: alternate sheet (η empty) --------------------------
Label(fig[2, 1], "(1b) η-plane: principal sheet only (v1)\n" *
                  "alt-sheet via ζ-frame above";
      tellwidth = false, halign = :center, valign = :center,
      fontsize = 11)

ax_ζ1 = Axis(fig[2, 2];
             title = "(2b) |w(ζ)|, ζ-plane sheet $alt_sheet_label",
             xlabel = "Re ζ", ylabel = "Im ζ",
             aspect = DataAspect(),
             limits = (ζ_RE_LO, ζ_RE_HI, ζ_IM_LO, ζ_IM_HI))
heatmap!(ax_ζ1, ζ_xs, ζ_ys, clip.(abs.(W_sheetalt));
         colormap = :viridis, colorrange = (0, ABS_CAP),
         nan_color = :transparent)
scatter!(ax_ζ1, [0.0], [0.0]; color = :orange, markersize = 10,
         marker = :circle, strokewidth = 1)
scatter!(ax_ζ1, [0.0], [2π]; color = :orange, markersize = 10,
         marker = :circle, strokewidth = 1)

ax_z1 = Axis(fig[2, 3];
             title = "(3b) log₁₀|u(z)|, z-plane sheet $alt_sheet_label (lift = $alt_lift)",
             xlabel = "Re z", ylabel = "Im z",
             aspect = DataAspect(),
             limits = (-Z_HALF, Z_HALF, -Z_HALF, Z_HALF))
heatmap!(ax_z1, z_xs, z_ys, Ulog_z_alt;
         colormap = :viridis, colorrange = (LOG_LO, LOG_HI),
         nan_color = :transparent)
lines!(ax_z1, [-Z_HALF, Z_HALF], [0.0, 0.0];
       color = (:black, 0.25), linewidth = 0.3)
lines!(ax_z1, [0.0, 0.0], [-Z_HALF, Z_HALF];
       color = (:black, 0.25), linewidth = 0.3)

Colorbar(fig[1, 4], hm_η; label = "|v| / |w|  (cap $ABS_CAP)")
Colorbar(fig[2, 4], hm_z; label = "log₁₀|u(z)|  ∈ [$LOG_LO, $LOG_HI]")

save(OUTPNG, fig)
@printf("\nWrote %s\n", OUTPNG)
