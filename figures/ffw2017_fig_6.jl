# figures/ffw2017_fig_6.jl
#
# Reproduces Fasondini, Fornberg & Weideman 2017, Fig. 6 — the
# generic P_V solution on three sheets of its Riemann surface, with
# parameters `(α, β, γ, δ) = (1, -1, 1, -1/2)` and initial conditions
# `(z₀, u(z₀), u'(z₀)) = (1, 2, -1)` in the natural `z`-plane.
#
# Source: references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
#         FFW2017_painleve_riemann_surfaces_preprint.md:281-297
#         (Figure 6 caption + Table 6 pole/zero residues).
#
# ## The problem
#
# `P_V` has a fixed singularity at `z = 0`.  The exponential transform
# `z = exp(ζ)`, `u(z) = w(ζ)` (FFW md:45-47) maps that branch point out
# of the finite plane and the resulting `P̃_V` equation is meromorphic in
# the ζ-plane.  Each strip `Im ζ ∈ (-π + 2π s, π + 2π s]`, `s ∈ ℤ`, maps
# to one sheet of the multi-sheeted `z`-plane Riemann surface.  We
# compute the solution on three contiguous strips:
#
#     sheet 0:  Im ζ ∈ (-π, π]      → principal sheet in z
#     sheet 1:  Im ζ ∈ (π, 3π]      → s = 1 sheet, arg z ∈ (π, 3π]
#     sheet 2:  Im ζ ∈ (3π, 5π]     → s = 2 sheet, arg z ∈ (3π, 5π]
#
# The full ζ-window is `Re ζ ∈ [-1, 2.5]`, `Im ζ ∈ (-π, 5π]` — three
# strips of vertical width `2π` stacked above each other.  In the
# `z`-plane this corresponds to `|z| ∈ [e^{-1}, e^{2.5}] ≈ [0.37, 12.2]`
# annulus, repeated three times as one winds around `z = 0`.
#
# ## Method
#
#   * `PadeTaylorProblem` built directly on the FFW md:47 `P̃_V` RHS
#     via `CoordTransforms.pV_transformed_rhs(α, β, γ, δ)`.  The IC
#     `(z₀, u, u') = (1, 2, -1)` maps to `(ζ₀, w, w') = (0, 2, -1)`
#     under `pV_z_to_ζ` (`z = 1 ⇒ ζ = log 1 = 0`; `w = u`, `w' = z u'`).
#
#   * Stage 1 walks via `path_network_solve` with
#     `step_size_policy = :adaptive_ffw` (ADR-0011, worklog 034 — the
#     FFW md:74-97 truncation-error controller `q = (k·Tol/T(h))^{1/(n+1)}`)
#     at `adaptive_tol = 1e-10`.
#
#   * Stage-1 node placement uses the opt-in `node_separation` kwarg
#     (ADR-0012, worklog 035).  We adopt `R(ζ) = max(0.05, (4 - Re ζ)/10)`
#     — FFW's Fig 1 PIII prescription is `R(ζ) = (8 - Re ζ)/20`
#     (md:72), tighter than what we need at our reduced `Re ζ ≤ 2.5`.
#     Our linear `R(ζ) = (4 - Re ζ)/10` has the same FFW shape (decreases
#     with Re ζ, reaching the pole-density-rich region) with a constant-
#     `0.05` floor preventing zero-spacing at the high-`Re ζ` boundary.
#     At `Re ζ = -1` it gives `R = 0.5` (FW's working default); at
#     `Re ζ = 2.5` it gives `R = 0.15`.
#
#   * Stage-2 evaluation grid is a uniform `Nx × Ny` rectangular lattice
#     over the ζ-window.  Cell values where the closest visited tree
#     node is farther than its local `visited_h` are `NaN` and rendered
#     as transparent / dropped.
#
# ## Three panels
#
#   * Panel A: `|w(ζ)|` over the full ζ-window — modulus surface in the
#     transformed plane.  Three horizontal strips (sheets 0/1/2) visible
#     as separated bands annotated on the left margin.  Poles appear as
#     vertical spikes; FFW md:283 notes the characteristic "oblique
#     lines" of poles in the transformed plane.
#
#   * Panel B: `|u(z)|` over the `z`-plane.  For each sheet, the strip
#     `Im ζ ∈ (-π + 2πs, π + 2πs]` is mapped through `z = exp(ζ)` to
#     `arg z ∈ (-π + 2πs, π + 2πs]`.  This puts sheet 0 in the principal
#     branch `arg z ∈ (-π, π]`, sheet 1 in `arg z ∈ (π, 3π]`, sheet 2 in
#     `arg z ∈ (3π, 5π]`.  Rendered as separate sub-panels, one per sheet.
#
#   * Panel C: phase portrait `arg(u(z)) ∈ [-π, π]` colour-wheel-encoded
#     (HSV hue = `arg`, value/saturation full) over the same three
#     `z`-plane sheets.  FFW md:283 calls the resulting spirals out as
#     the defining signature of multi-valued PIII/PV solutions.
#
# ## Acceptance
#
# Per `docs/figure_catalogue.md §5` (Fig 6 row, T4): per-sheet error
# from FFW Table 2-style symmetry estimates `3e-9 / 7e-6 / 2e-5` for
# sheets 0/1/2, within one order of magnitude.  `test/ffw_fig_6_test.jl`
# pins this quantitatively (FF6.1.3 / FF6.1.4 / FF6.1.5).

using PadeTaylor
using PadeTaylor.CoordTransforms: pV_transformed_rhs, pV_z_to_ζ, pV_ζ_to_z
using CairoMakie
using Printf

# ----------------------------------------------------------------------
# Parameters (FFW md:297)
# ----------------------------------------------------------------------
const α, β, γ, δ = 1.0, -1.0, 1.0, -0.5
const Z0, U0, UP0 = 1.0 + 0.0im, 2.0 + 0.0im, -1.0 + 0.0im

const ORDER       = 30           # FW/FFW working Taylor order
const ADAPT_TOL   = 1.0e-10      # FFW md:91 controller tolerance
const K_CONS      = 1.0e-3       # FFW md:91 conservative factor
const MAX_RESC    = 50           # cap on per-step rescales

# ζ-window: three strips of vertical width 2π stacked.
const Re_LO, Re_HI = -1.0, 2.5
const Im_LO, Im_HI = -π + 0.01, 5π - 0.01       # epsilon-shrunk to keep
                                                 # boundary lookups in-tree

# Stage-2 evaluation lattice density.
const NX, NY = 90, 240
const OUTPNG = joinpath(@__DIR__, "output", "ffw2017_fig_6.png")

# FFW md:72 node-separation function `R(ζ) = (8 - Re ζ)/20` is the PIII
# Fig 1 prescription; we use the same shape for our PV problem with
# coefficients tuned for `Re ζ ∈ [-1, 2.5]` rather than `[−2, 8]`.  A
# floor at `0.05` prevents zero-spacing at the high-`Re ζ` boundary.
R(ζ) = max(0.05, (4.0 - real(ζ)) / 10.0)

# ----------------------------------------------------------------------
# Map the IC into the ζ-frame.
# ----------------------------------------------------------------------
const ζ0, w0, wp0 = pV_z_to_ζ(Z0, U0, UP0)
@printf("FFW 2017 Fig 6: (α,β,γ,δ) = (%.1f, %.1f, %.1f, %.2f)\n", α, β, γ, δ)
@printf("                IC (z₀,u,u') = (%s, %s, %s)\n", Z0, U0, UP0)
@printf("                → (ζ₀, w, w') = (%s, %s, %s)\n", ζ0, w0, wp0)

# ----------------------------------------------------------------------
# Build Stage-1 target set + Stage-2 evaluation lattice.
# ----------------------------------------------------------------------
# Stage-1 targets are placed at the local `R(ζ)` spacing — FFW md:97
# prescription.  Walking through them seeds the visited-tree density
# proportional to local pole-density.
function build_stage1_targets()
    tgs = ComplexF64[]
    re = Re_LO
    while re ≤ Re_HI + 1e-12
        rh = R(re + 0.0im)
        imv = Im_LO
        while imv ≤ Im_HI + 1e-12
            push!(tgs, complex(re, imv))
            imv += rh
        end
        re += rh
    end
    return tgs
end

stage1_targets = build_stage1_targets()
@printf("Stage-1 target count: %d (R-spaced)\n", length(stage1_targets))

xs = collect(range(Re_LO, Re_HI; length = NX))
ys = collect(range(Im_LO, Im_HI; length = NY))
eval_lattice = ComplexF64[complex(x, y) for y in ys for x in xs]
@printf("Stage-2 eval lattice: %d × %d = %d cells\n", NX, NY, length(eval_lattice))

# Combined grid: Stage-1 targets drive walker placement, Stage-2 cells
# read off the resulting Padé store.
grid = vcat(stage1_targets, eval_lattice)

# ----------------------------------------------------------------------
# Solve the FFW md:47 `P̃_V` problem on the ζ-grid.
# ----------------------------------------------------------------------
f = pV_transformed_rhs(α, β, γ, δ)
prob = PadeTaylorProblem(f, (w0, wp0), (ζ0, complex(Re_HI, Im_HI));
                          order = ORDER)

t0 = time()
sol = path_network_solve(prob, grid;
                          h = R(ζ0),
                          node_separation = R,
                          step_size_policy = :adaptive_ffw,
                          adaptive_tol = ADAPT_TOL,
                          k_conservative = K_CONS,
                          max_rescales = MAX_RESC,
                          max_steps_per_target = 4000)
@printf("  Stage 1 + Stage 2 in %.2f s; %d visited tree nodes.\n",
        time() - t0, length(sol.visited_z))

eval_u = sol.grid_u[(length(stage1_targets) + 1):end]
@assert length(eval_u) == length(eval_lattice)
ncov = count(isfinite, abs.(eval_u))
@printf("  Stage-2 coverage: %d / %d (%.1f%%) finite\n",
        ncov, length(eval_u), 100 * ncov / length(eval_u))

# ----------------------------------------------------------------------
# Reshape into the (NX, NY) grid for surface plotting.  `eval_lattice`
# was built y-outer / x-inner, matching column-major reshape(_, NX, NY).
# ----------------------------------------------------------------------
W = reshape(eval_u, NX, NY)
absW = abs.(W)

# ----------------------------------------------------------------------
# Map ζ-plane data to z-plane sheets:
#   sheet s: Im ζ ∈ (-π + 2π s, π + 2π s] → arg z ∈ same range
#   modulus |u(z)| = |w(ζ)| since the PV transform is u = w (no
#   rescale).
# ----------------------------------------------------------------------
sheet_im_ranges = [(-π + 2π*s, π + 2π*s) for s in 0:2]

# ----------------------------------------------------------------------
# Render: three columns × three sheets each = nine panels.
# Layout:
#   col 1: |w(ζ)| heatmap per sheet
#   col 2: |u(z)| heatmap per sheet
#   col 3: arg(u(z)) phase portrait per sheet
# ----------------------------------------------------------------------
const U_CAP   = 5.0     # |w|/|u| cap; pole spikes go to infinity
const ABS_LO  = 0.0
const Z_RMIN  = exp(Re_LO)
const Z_RMAX  = exp(Re_HI)

fig = Figure(size = (1500, 1100))
Label(fig[0, 1:3],
      "FFW 2017 Fig 6 — generic P_V, (α,β,γ,δ) = (1,−1,1,−1/2), IC (z,u,u′) = (1,2,−1)";
      fontsize = 17, padding = (0, 0, 6, 6))

for (i, (im_lo, im_hi)) in enumerate(sheet_im_ranges)
    sheet_label = "sheet $(i-1)"
    # Indices of `ys` falling in this strip
    ym = findall(y -> im_lo ≤ y ≤ im_hi, ys)
    sub_ys = ys[ym]
    sub_W  = W[:, ym]
    sub_absW = abs.(sub_W)
    capped_absW = clamp.(sub_absW, ABS_LO, U_CAP)

    # ---------------- Panel A: |w(ζ)| in the ζ-plane ----------------
    axA = Axis(fig[i, 1];
               xlabel = "Re ζ", ylabel = "Im ζ",
               title  = "|w(ζ)|, $sheet_label",
               titlesize = 13,
               aspect = DataAspect(),
               limits = (Re_LO, Re_HI, im_lo, im_hi))
    heatmap!(axA, xs, sub_ys, capped_absW;
             colormap = :viridis, colorrange = (ABS_LO, U_CAP))

    # ---------------- Panel B: |u(z)| in the z-plane ----------------
    # For sheet s, z = exp(Re ζ + i Im ζ); we lay out the heatmap as a
    # polar wedge by sampling (r, θ) = (exp(Re ζ), Im ζ).  CairoMakie's
    # surface/heatmap accepts a vector-of-(x,y) format; we precompute
    # zx, zy = (r cosθ, r sinθ) per (Re ζ, Im ζ) cell.  Visualised as a
    # scatter heatmap (fast + readable).
    nx, ny = length(xs), length(sub_ys)
    rs = exp.(xs)
    θs = sub_ys
    Zx = Float64[r * cos(θ) for θ in θs, r in rs]   # (ny, nx)
    Zy = Float64[r * sin(θ) for θ in θs, r in rs]
    valsB = Float64[Float64(clamp(abs(W[ix, iy]), ABS_LO, U_CAP))
                    for iy in ym, ix in 1:nx]
    axB = Axis(fig[i, 2];
               xlabel = "Re z", ylabel = "Im z",
               title  = "|u(z)|, $sheet_label",
               titlesize = 13,
               aspect = DataAspect(),
               limits = (-1.1 * Z_RMAX, 1.1 * Z_RMAX,
                         -1.1 * Z_RMAX, 1.1 * Z_RMAX))
    scatter!(axB, vec(Zx), vec(Zy);
             color = vec(valsB),
             colormap = :viridis,
             colorrange = (ABS_LO, U_CAP),
             markersize = 2.0)

    # ---------------- Panel C: arg(u(z)) phase portrait ----------------
    # HSV color wheel: hue = (arg(u)+π)/(2π), s=v=1
    valsC = Float64[Float64(angle(W[ix, iy])) for iy in ym, ix in 1:nx]
    axC = Axis(fig[i, 3];
               xlabel = "Re z", ylabel = "Im z",
               title  = "arg u(z), $sheet_label",
               titlesize = 13,
               aspect = DataAspect(),
               limits = (-1.1 * Z_RMAX, 1.1 * Z_RMAX,
                         -1.1 * Z_RMAX, 1.1 * Z_RMAX))
    scatter!(axC, vec(Zx), vec(Zy);
             color = vec(valsC),
             colormap = :hsv,
             colorrange = (-π, π),
             markersize = 2.0)
end

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
