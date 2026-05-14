# figures/fw2011_fig_5_2.jl
#
# Reproduces Fornberg & Weideman 2011, Fig. 5.2 — the accuracy/cost
# trade-off of the Padé integrator as a function of the two free
# parameters `order` and `h`, which is what motivates FW 2011's
# standard choice `order = 30`, `h = 0.5`.
#
# Source: references/markdown/FW2011_painleve_methodology_JCP230/
#         FW2011_painleve_methodology_JCP230.md:277-326.
#
# The test problem is the Weierstrass ℘ medium-range integration of
# `fw2011_fig_5_1.jl`: u'' = 6 u^2 with the equianharmonic IC
#
#     u(0)  = 1.071822516416917,
#     u'(0) = 1.710337353176786,
#
# integrated to z = 30, whose exact value u(30) = 1.095098255959744
# (FW eq. 5.3) gives a clean relative-error oracle.
#
# Two panels (FW md:320-326):
#
#   (a) log10(relative error) as a surface over the (order, h) plane.
#       FW truncate the surface where the relative error exceeds a cap
#       (their figure caption value is clipped in the markdown; we cap
#       log10(rel err) at 0, i.e. rel err = 1, which is the natural
#       "no useful accuracy" plateau and also catches the (large-h,
#       low-order) corner where the integrator cannot bridge a pole).
#
#   (b) contour lines of equal accuracy (solid) and of equal compute
#       time (dashed).  Their intersections single out the favourable
#       region around order = 30, h = 0.5.
#
# Each (order, h) cell is one `path_network_solve` walk from z = 0 to
# z = 30 — the same routine `test/pathnetwork_test.jl` exercises under
# PN.2.2.  A walk that throws (the integrator genuinely failing for
# that parameter pair) is recorded as rel err = 1, which is exactly
# the truncation plateau.

using PadeTaylor
using CairoMakie
using Printf

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
const Z_END   = 30.0                       # medium-range target (md:299)
const U0      = 1.071822516416917          # ℘ IC, FW md:294
const UP0     = 1.710337353176786          # ℘ IC, FW md:295
const U30_REF = 1.095098255959744          # exact u(30), FW eq. 5.3
const ORDERS  = 4:2:50                     # FW Fig 5.2 order axis
const HS      = 0.05:0.05:1.0              # FW Fig 5.2 h axis
const ERR_CAP = 0.0                        # log10(rel err) truncation cap
const ERR_FLOOR = -16.0                    # log10(rel err) display floor
const OUTPNG  = joinpath(@__DIR__, "output", "fw2011_fig_5_2.png")

# Weierstrass ℘ test problem: u'' = 6 u^2 (FW eq. 5.1).
pweier(z, u, up) = 6 * u^2

# ----------------------------------------------------------------------
# Sweep the (order, h) plane
# ----------------------------------------------------------------------
nO, nH = length(ORDERS), length(HS)
logerr = fill(NaN, nO, nH)        # logerr[i, j] at (ORDERS[i], HS[j])
walltime = fill(NaN, nO, nH)

@printf("FW 2011 Fig 5.2: sweeping %d orders × %d h-values = %d ℘ walks to z=30 ...\n",
        nO, nH, nO * nH)
t_sweep = time()
for (j, h) in enumerate(HS), (i, ord) in enumerate(ORDERS)
    prob = PadeTaylorProblem(pweier, (U0, UP0), (0.0, Z_END); order = ord)
    maxsteps = max(400, round(Int, 4 * Z_END / h))
    t0 = time()
    relerr = try
        # NB: `path_network_solve` carries its OWN `order` kwarg and
        # ignores `prob.order` — the sweep variable must be passed
        # here, not (only) to `PadeTaylorProblem`.
        sol = path_network_solve(prob, ComplexF64[Z_END + 0im];
                                 h = h, order = ord,
                                 max_steps_per_target = maxsteps)
        abs(sol.grid_u[1] - U30_REF) / abs(U30_REF)
    catch
        # The integrator genuinely failed for this (order, h): record
        # it at the truncation plateau (rel err = 1).
        1.0
    end
    walltime[i, j] = time() - t0
    logerr[i, j]   = clamp(log10(max(relerr, 1e-300)), ERR_FLOOR, ERR_CAP)
end
@printf("  swept in %.1f s; best log10(rel err) = %.2f at (order, h) = ",
        time() - t_sweep, minimum(logerr))
let ibest = argmin(logerr)
    @printf("(%d, %.2f)\n", ORDERS[ibest[1]], HS[ibest[2]])
end
# Accuracy + time at FW's standard choice (order = 30, h = 0.5).
let i30 = findfirst(==(30), ORDERS), j05 = findfirst(==(0.5), HS)
    @printf("  at FW's standard (order=30, h=0.5): log10(rel err) = %.2f, wall = %.3f s\n",
            logerr[i30, j05], walltime[i30, j05])
end

orders = collect(Float64, ORDERS)
hs     = collect(Float64, HS)

# 3×3 box smoother — FW md:326 explicitly smooths the raw sweep data
# before drawing the Fig 5.2b contours ("some fine-scale 'jitter' …
# was smoothed out when we created these contours").  Panel (a) keeps
# the raw surface; panel (b) contours the smoothed fields.
function smooth3(M)
    nr, nc = size(M)
    S = similar(M)
    for j in 1:nc, i in 1:nr
        acc, cnt = 0.0, 0
        for dj in -1:1, di in -1:1
            ii, jj = i + di, j + dj
            if 1 ≤ ii ≤ nr && 1 ≤ jj ≤ nc && isfinite(M[ii, jj])
                acc += M[ii, jj]; cnt += 1
            end
        end
        S[i, j] = cnt == 0 ? M[i, j] : acc / cnt
    end
    return S
end

# ----------------------------------------------------------------------
# Render: (a) the error surface, (b) accuracy + time contours
# ----------------------------------------------------------------------
fig = Figure(size = (1180, 470))

# Panel (a): log10(rel err) surface over the (order, h) plane.
ax_a = Axis3(fig[1, 1];
             xlabel = "order", ylabel = "h", zlabel = "log₁₀(relative error)",
             title  = "(a)  accuracy vs order and h",
             azimuth = 1.30π, elevation = 0.18π, aspect = (1.0, 1.0, 0.6))
surface!(ax_a, orders, hs, logerr; colormap = :viridis)
wireframe!(ax_a, orders, hs, logerr; color = (:black, 0.25), linewidth = 0.4)

# Panel (b): equal-accuracy contours (solid) and equal-time contours
# (dashed); their intersection picks out the favourable region.
# Both fields are 3×3-smoothed first, per FW md:326.
ax_b = Axis(fig[1, 2];
            xlabel = "order", ylabel = "h",
            title  = "(b)  accuracy (solid) and compute-time (dashed) contours",
            limits = (first(orders), last(orders), first(hs), last(hs)))
contour!(ax_b, orders, hs, smooth3(logerr);
         levels = Float64[-2, -4, -6, -8, -10, -12],
         color = :black, linewidth = 1.6, labels = true, labelsize = 10)
# Compute time, normalised to its minimum over the swept plane.
reltime = walltime ./ minimum(filter(isfinite, walltime))
contour!(ax_b, orders, hs, smooth3(reltime);
         levels = Float64[2, 5, 15, 50], color = :gray45,
         linewidth = 1.2, linestyle = :dash)
# FW's standard choice.
scatter!(ax_b, [30.0], [0.5]; color = :red, markersize = 13)
text!(ax_b, 30.5, 0.52; text = "order=30, h=0.5", color = :red, fontsize = 11)

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
