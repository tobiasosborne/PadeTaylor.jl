# figures/fw2011_fig_3_3.jl
#
# Reproduces Fornberg & Weideman 2011, Fig. 3.3 — the pole-field edge
# detector: log10|Δu| obtained by applying the 5-point Laplacian
# stencil (FW 2011 eq. 3.3) to a numerical near-tritronquée PI
# solution, with the level-0.001 contour overlaid as the pole-field
# edge description.
#
# Source: references/markdown/FW2011_painleve_methodology_JCP230/
#         FW2011_painleve_methodology_JCP230.md:202-208.
#
# The idea (FW 2011 Section 3.2.2).  Away from its poles a meromorphic
# u(z) = u(x, y) is harmonic: u_xx + u_yy = 0.  The residual of the
# 5-point cross stencil
#
#     Δu ≈ ( u(x-h,y) + u(x+h,y) + u(x,y-h) + u(x,y+h) - 4 u(x,y) ) / h²
#
# is therefore O(h²) and tiny in smooth regions, but large inside pole
# fields where u is far from harmonic.  Plotting log10|Δu| gives a
# surface that is a deep flat plain over the smooth region and a range
# of sharp ridges over the pole fields; a single contour level cleanly
# separates the two.  FW uses level 0.001 (in the plotted log10|Δu|
# units — i.e. |Δu| ≈ 1.002), which is the default of
# `PadeTaylor.EdgeDetector.pole_field_mask`.
#
# This figure is the visual showcase of the `EdgeDetector` module
# (Phase 12.5): the stencil itself is `EdgeDetector.laplacian_residual`
# verbatim, and the contour level is `pole_field_mask`'s default.
#
# Method.  FW md:198 captions Fig. 3.3 "applied to a numerical
# near-tritronquee solution in the vicinity of the origin (cf. Figs.
# 1.1c and 3.1)".  Fig. 1.1c is the *tritronquee* case — the limiting
# solution with four of the five sectors entirely pole-free — and
# FW's Fig. 3.3 shows exactly that signature: one pole field in the
# back-right of the window and a single large clean smooth plain
# elsewhere.  We therefore use the precise tritronquee initial
# conditions FW later pin in eq. (4.1) (md:224-227):
#
#     u(0)  = -0.1875543083404949,
#     u'(0) =  0.3049055602612289.
#
# (The cruder Fig.-3.1 ICs -0.1875 / 0.3049 instead place displaced
# pole fields in four sectors, which would fill this window with far
# more pole-field area than FW's figure shows.)  The solution is
# sampled on a uniform lattice over the near-origin window
# x ∈ [-4, 8], y ∈ [-6, 6] with `path_network_solve` and Schwarz
# symmetry (PI has real coefficients, the ICs are on the real axis),
# then handed to `laplacian_residual`.
#
# FW 2011 (md:208) explicitly notes that "low level 'ridges' appear
# within the flat areas" because neighbouring lattice entries can be
# reached by entirely different paths through the complex plane, then
# amplified by the 1/h² in the stencil.  Expect — and do not be
# alarmed by — a few such ridges in the smooth plain here too.

using PadeTaylor
using PadeTaylor.EdgeDetector: laplacian_residual, pole_field_mask
using CairoMakie
using Printf

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
const X_LO, X_HI = -4.0, 8.0     # near-origin window, FW Fig 3.3 axes
const Y_LO, Y_HI = -6.0, 6.0
const DGRID  = 0.25              # lattice spacing fed to the 5-pt stencil
const H_STEP = 0.5               # FW 2011 path-network step length (md:164)
const ORDER  = 30                # FW 2011 Section 5.1 Taylor order
const LEVEL  = 0.001             # FW Fig 3.3 contour level (log10|Δu| units)
const Z_LO, Z_HI = -7.0, 5.0     # log10|Δu| display clamp, matching FW framing
const OUTPNG = joinpath(@__DIR__, "output", "fw2011_fig_3_3.png")

# Precise tritronquée initial conditions, FW 2011 eq. (4.1) (md:226).
const U0  = -0.1875543083404949
const UP0 =  0.3049055602612289

# PI right-hand side: u'' = 6 u^2 + z.
pI(z, u, up) = 6 * u^2 + z

# ----------------------------------------------------------------------
# Solve PI on the lattice
# ----------------------------------------------------------------------
xs = range(X_LO, X_HI; step = DGRID)
ys = range(Y_LO, Y_HI; step = DGRID)
nx, ny = length(xs), length(ys)
# Flatten x-major so that `reshape(_, ny, nx)` yields u_grid[iy, ix] —
# the rows-index-y / cols-index-x orientation `laplacian_residual`
# documents (EdgeDetector module docstring).
grid = ComplexF64[complex(x, y) for x in xs for y in ys]

prob = PadeTaylorProblem(pI, (U0, UP0), (0.0, 10.0); order = ORDER)

@printf("FW 2011 Fig 3.3: solving PI on a %d x %d lattice, x∈[%g,%g] y∈[%g,%g] ...\n",
        nx, ny, X_LO, X_HI, Y_LO, Y_HI)
t0  = time()
sol = path_network_solve(prob, grid;
                         h = H_STEP,
                         max_steps_per_target = 4000,
                         enforce_real_axis_symmetry = true)
u_grid = reshape(sol.grid_u, ny, nx)          # u_grid[iy, ix]
nbad   = count(!isfinite, u_grid)
@printf("  done in %.1f s; %d visited tree nodes; %d / %d lattice points finite.\n",
        time() - t0, length(sol.visited_z), nx * ny - nbad, nx * ny)

# ----------------------------------------------------------------------
# 5-point Laplacian stencil (FW 2011 eq. 3.3) via the EdgeDetector module
# ----------------------------------------------------------------------
Δu   = laplacian_residual(u_grid, DGRID)      # boundary cells -> NaN
mask = pole_field_mask(u_grid, DGRID; level = LEVEL)
L    = log10.(abs.(Δu))                       # L[iy, ix]; boundary NaN

# How cleanly does the level-0.001 contour separate the two regimes?
finite_L = filter(isfinite, vec(L))
@printf("  log10|Δu| over interior: min %.2f, max %.2f, median %.2f\n",
        minimum(finite_L), maximum(finite_L),
        sort(finite_L)[fld(length(finite_L), 2)])
@printf("  pole_field_mask (level = %.3f): %d / %d interior cells flagged as pole-field.\n",
        LEVEL, count(mask), (nx - 2) * (ny - 2))

# ----------------------------------------------------------------------
# Render: log10|Δu| surface (interior only) + the level-LEVEL contour,
# drawn on the surface as FW does.
# ----------------------------------------------------------------------
# `surface` / `contour` want Z[ix, iy]; L is [iy, ix].  Trim the 1-cell
# NaN border the stencil leaves, then transpose.
xi = xs[2:end-1]
yi = ys[2:end-1]
Z  = clamp.(permutedims(L[2:end-1, 2:end-1]), Z_LO, Z_HI)   # Z[ix, iy]

fig = Figure(size = (1000, 660))
ax  = Axis3(fig[1, 1];
            xlabel = "x", ylabel = "y", zlabel = "log₁₀|Δu|",
            title  = "FW 2011 Fig. 3.3 — pole-field edge detector, level-$(LEVEL) contour",
            azimuth = 1.30π, elevation = 0.16π,
            aspect = (1.0, 1.0, 0.5),
            zticks = -6:2:4)
surface!(ax, xi, yi, Z; colormap = :viridis, colorrange = (Z_LO, Z_HI))
wireframe!(ax, xi, yi, Z; color = (:black, 0.18), linewidth = 0.3)
# The pole-field edge: the level-LEVEL contour of log₁₀|Δu|, drawn at
# its true height on the surface (FW 2011's "level curve 0.001").
contour!(ax, xi, yi, Z; levels = [Float64(LEVEL)],
         color = :black, linewidth = 2.5)

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
