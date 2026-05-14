# figures/fw2011_fig_4_2.jl
#
# Reproduces Fornberg & Weideman 2011, Fig. 4.2 — PI solutions along
# the real axis obeying u(0) = 0.
#
#     u''(z) = 6 u(z)^2 + z
#
# Two panels (FW md:239-241):
#
#   (a) one tronquée case bracketed by two near-tronquée solutions:
#       u'(0) = 1.8518            (just below the tronquée value)
#       u'(0) = 1.85185403375822  (the tronquée case)
#       u'(0) = 1.8519            (just above the tronquée value)
#   (b) the two tronquée cases obeying u(0) = 0:
#       u'(0) =  1.85185403375822
#       u'(0) = -0.45142740474177
#
# In both panels the dotted curves are the leading-term asymptotics
# ±√(-z/6) of FW eq. (1.2), valid for z = x < 0.
#
# Source: references/markdown/FW2011_painleve_methodology_JCP230/
#         FW2011_painleve_methodology_JCP230.md:231-241.
#
# The point of the figure (FW md:233): as u'(0) crosses the tronquée
# value 1.85185403375822, the left pole field is pushed out along the
# negative real axis and returns with its front row of poles flipping
# between straddling the axis and sitting on it.  The near-tronquée
# curves (1.8518, 1.8519) therefore carry real-axis poles — visible as
# the near-vertical asymptotes — while the tronquée curve (panel b
# solid) stays smooth and hugs the lower leading-term branch.
#
# Method.  FW md:233 notes these real-axis curves ARE the thick lines
# drawn on the Fig. 4.3 / 4.4 pole-field surfaces — i.e. FW reads them
# off the 2-D path network rather than stepping straight along the
# real axis through the poles.  We do the same: `path_network_solve`
# on a 1-D grid of real-axis target points.  The wedge walk detours
# into the complex plane around each pole, and Stage 2 evaluates the
# local Padé back down at the real target — large but finite where the
# target sits on a pole, which is exactly the asymptote we want.

using PadeTaylor
using CairoMakie
using Printf

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
const H_STEP = 0.5               # FW 2011 working step length (md:164)
const ORDER  = 30                # FW 2011 Section 5.1 Taylor order
const XS     = range(-10.5, 1.5; length = 241)   # real-axis sample grid
const U_CLAMP = 50.0             # clip |u| before plotting (poles -> ±∞)
const OUTPNG = joinpath(@__DIR__, "output", "fw2011_fig_4_2.png")

# The tronquée and near-tronquée initial slopes, FW Fig 4.2 legends.
const UP0_TRO_A = 1.85185403375822    # tronquée  (panels a + b, solid)
const UP0_BELOW = 1.8518              # near-tronquée, just below (a, dashed)
const UP0_ABOVE = 1.8519              # near-tronquée, just above (a, dash-dot)
const UP0_TRO_B = -0.45142740474177   # the other tronquée case (b, dashed)

# PI right-hand side: u'' = 6 u^2 + z.
pI(z, u, up) = 6 * u^2 + z

# ----------------------------------------------------------------------
# Real-axis profile u(x) for a given initial slope, via the path network
# ----------------------------------------------------------------------
"""
    realaxis_u(up0) -> Vector{Float64}

`u(x)` sampled on `XS` for the PI solution with `u(0) = 0`,
`u'(0) = up0`.  `path_network_solve` covers the real-axis target grid
by detouring around the poles; the returned values are `real(u)`
clipped to `±U_CLAMP` so the renderer is not handed `Inf` at a pole.
"""
function realaxis_u(up0::Real)
    grid = ComplexF64[complex(x, 0.0) for x in XS]
    prob = PadeTaylorProblem(pI, (0.0, float(up0)), (0.0, 10.0); order = ORDER)
    t0   = time()
    sol  = path_network_solve(prob, grid;
                              h = H_STEP, max_steps_per_target = 6000)
    u    = clamp.(real.(sol.grid_u), -U_CLAMP, U_CLAMP)
    npole = count(>(5.0), abs.(u))
    @printf("  u'(0) = %-18.14g : %.1f s, %d visited nodes, %d near-pole samples\n",
            up0, time() - t0, length(sol.visited_z), npole)
    return u
end

# Leading-term asymptotics ±√(-x/6) of FW eq. (1.2); real only for x<0.
lead_plus(x)  = x < 0 ? sqrt(-x / 6) : NaN
lead_minus(x) = x < 0 ? -sqrt(-x / 6) : NaN

# ----------------------------------------------------------------------
# Compute every curve
# ----------------------------------------------------------------------
println("FW 2011 Fig 4.2: integrating PI real-axis profiles ...")
u_below = realaxis_u(UP0_BELOW)
u_tro_a = realaxis_u(UP0_TRO_A)
u_above = realaxis_u(UP0_ABOVE)
u_tro_b = realaxis_u(UP0_TRO_B)

xs   = collect(XS)
lp   = lead_plus.(xs)
lm   = lead_minus.(xs)

# ----------------------------------------------------------------------
# Render: two panels, FW's [-10, 2] × [-3, 2] framing
# ----------------------------------------------------------------------
fig = Figure(size = (1100, 470))

function panel(pos, title)
    ax = Axis(fig[1, pos]; xlabel = "x", ylabel = "u(x)", title = title,
              limits = (-10, 2, -3, 2))
    # Leading-term ±√(-x/6) reference branches.
    lines!(ax, xs, lp; color = :gray35, linestyle = :dot, linewidth = 1.2)
    lines!(ax, xs, lm; color = :gray35, linestyle = :dot, linewidth = 1.2)
    # FW's dotted axes through the origin.
    hlines!(ax, [0.0]; color = :gray70, linestyle = :dot, linewidth = 0.8)
    vlines!(ax, [0.0]; color = :gray70, linestyle = :dot, linewidth = 0.8)
    return ax
end

# Panel (a): tronquée bracketed by the two near-tronquée solutions.
ax_a = panel(1, "(a)")
l_below = lines!(ax_a, xs, u_below; color = :black, linestyle = :dash,    linewidth = 1.5)
l_tro_a = lines!(ax_a, xs, u_tro_a; color = :black, linestyle = :solid,   linewidth = 1.5)
l_above = lines!(ax_a, xs, u_above; color = :black, linestyle = :dashdot, linewidth = 1.5)
axislegend(ax_a,
    [l_below, l_tro_a, l_above],
    ["u'(0) = 1.8518", "u'(0) = 1.85185403375822", "u'(0) = 1.8519"];
    position = :rb, framevisible = true, labelsize = 10)

# Panel (b): the two tronquée cases obeying u(0) = 0.
ax_b = panel(2, "(b)")
l_tro_b1 = lines!(ax_b, xs, u_tro_a; color = :black, linestyle = :solid, linewidth = 1.5)
l_tro_b2 = lines!(ax_b, xs, u_tro_b; color = :black, linestyle = :dash,  linewidth = 1.5)
axislegend(ax_b,
    [l_tro_b1, l_tro_b2],
    ["u'(0) =  1.85185403375822", "u'(0) = -0.45142740474177"];
    position = :rb, framevisible = true, labelsize = 10)

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
