# figures/fw2011_fig_4_8.jl
#
# Reproduces Fornberg & Weideman 2011, Fig. 4.8 — the Painlevé-I
# pole field for the initial conditions
#
#     u(0) = -5,    u'(0) = 0,
#
# displayed over the rectangle z = x + iy, -90 <= x <= 30,
# -30 <= y <= 30.
#
# Source: references/markdown/FW2011_painleve_methodology_JCP230/
#         FW2011_painleve_methodology_JCP230.md:265-314.
#
# FW md:271 (§4.5): "Fig. 4.8 illustrates that changes in pole patterns
# need neither be gradual, nor be associated with well defined smooth
# bands.  With the ICs u(0) = -5, u'(0) = 0, a sharp pattern transition
# is visible near Re(z) = -60."  The large IC produces, near the
# origin, a six-fold rather than five-fold local pattern (FW md:265);
# moving out along the negative real axis the pole lattice abruptly
# re-organises into the asymptotic five-fold structure.  The reproduction
# target (figure_catalogue §1): the transition's Re-coordinate matches
# FW's ≈ -60 to within ±0.2 (visual).
#
# ----------------------------------------------------------------------
# Method
#
# Same as `fw2011_fig_4_7.jl`: `edge_gated_pole_field_solve` (the
# region-growing IVP confined to the pole field, FW md:401, bead
# `padetaylor-dmb`) followed by `PoleField.extract_poles` (bead
# `padetaylor-xvf`).  FW computed Fig. 4.8 as "a composite of 18
# independent runs" (md:147) — a pure IVP fill, because this solution
# carries poles throughout the window with no extended smooth bands.
# The edge-gated solve reduces to exactly that pure IVP fill when the
# pole field covers the window, while still confining the path-network
# safely around the sharp transition region.
#
# Step length h = 0.5, Taylor order = 30 — FW 2011 working defaults
# (md:279).

using PadeTaylor
using CairoMakie
using Printf
include(joinpath(@__DIR__, "figutil.jl"))

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
const XLIM   = (-90.0, 30.0)     # FW md:314 window, real axis
const YLIM   = (-30.0, 30.0)     # FW md:314 window, imaginary axis
const STEP   = 1.0               # uniform lattice step — fine enough for the
                                 # §3.2.2 edge detector (EdgeGatedSolve
                                 # "Grid resolution matters")
const H_STEP = 0.5               # FW 2011 working step length (md:279)
const ORDER  = 30                # FW 2011 standard Taylor order (md:279)
const GROW   = 4                 # edge-gated region-growing ring width
const U0     = -5.0              # FW Fig 4.8 IC (md:314)
const UP0    =  0.0              # FW Fig 4.8 IC (md:314)
const X_TRANS = -60.0            # FW's stated transition location (md:271)
const OUTPNG = joinpath(@__DIR__, "output", "fw2011_fig_4_8.png")

# PI right-hand side: u'' = 6 u^2 + z.
pI(z, u, up) = 6 * u^2 + z

# ----------------------------------------------------------------------
# Solve + extract poles
# ----------------------------------------------------------------------
xs = range(XLIM[1], XLIM[2]; step = STEP)
ys = range(YLIM[1], YLIM[2]; step = STEP)
@printf("FW 2011 Fig 4.8 — PI pole field, u(0)=%g, u'(0)=%g, over [%g,%g]×[%g,%g]\n",
        U0, UP0, XLIM..., YLIM...)
@printf("  lattice %d × %d, step %g\n", length(xs), length(ys), STEP)

prob = PadeTaylorProblem(pI, (U0, UP0), (0.0, abs(XLIM[1])); order = ORDER)

t0  = time()
egs = edge_gated_pole_field_solve(prob, xs, ys;
                                  h = H_STEP, order = ORDER, grow_rings = GROW)
poles = extract_poles(egs.pn_solution)
@printf("  %d region-growing passes, %d pole-field cells, %d visited nodes; %d poles extracted; %.1f s\n",
        egs.iterations, count(egs.field_mask),
        length(egs.pn_solution.visited_z), length(poles), time() - t0)

# Quantitative companion to the visual acceptance: report the gap in the
# pole field along the negative real axis.  FW md:271 places the sharp
# transition near Re(z) = -60; we print the |z| range of the on-axis
# poles bracketing X_TRANS so the reproduction can be checked.
on_axis = sort!(Float64[real(p) for p in poles if abs(imag(p)) < STEP/2])
left_of  = isempty(filter(<(X_TRANS), on_axis)) ? NaN : maximum(filter(<(X_TRANS), on_axis))
right_of = isempty(filter(>(X_TRANS), on_axis)) ? NaN : minimum(filter(x -> X_TRANS < x < 0, on_axis))
@printf("  on-axis poles bracketing Re=%.0f: nearest below = %.2f, nearest above = %.2f\n",
        X_TRANS, left_of, right_of)

# ----------------------------------------------------------------------
# Render
# ----------------------------------------------------------------------
fig = Figure(; size = (1180, 640))
ax  = pole_scatter_axis(fig[1, 1], poles;
                        title = "FW 2011 Fig. 4.8 — PI pole field, u(0) = -5, u'(0) = 0",
                        xlim = XLIM, ylim = YLIM, markersize = 2.2)
# FW's stated transition location, md:271.
vlines!(ax, [X_TRANS]; color = (:red, 0.5), linewidth = 1.0, linestyle = :dash)
text!(ax, X_TRANS, YLIM[2] - 3; text = "Re(z) ≈ -60  (FW md:271)",
      color = (:red, 0.7), fontsize = 11, align = (:center, :top))

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
