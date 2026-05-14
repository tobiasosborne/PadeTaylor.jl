# figures/fw2011_fig_5_1.jl
#
# Reproduces Fornberg & Weideman 2011, Fig. 5.1 — the Weierstrass-℘
# test problem: the poles of the ℘-function (large dots) together with
# the path the Padé integrator takes to reach z = 30.
#
# Source: references/markdown/FW2011_painleve_methodology_JCP230/
#         FW2011_painleve_methodology_JCP230.md:281-318.
#
# FW 2011 §5.1.1 drops the z-term of PI to get a test problem with the
# same in-pole-field numerics but a *known* solution:
#
#     u''(z) = 6 u(z)^2,        u(z) = ℘(z + c₁; 0, c₂),
#
# with c₁ = -1, c₂ = 2, which fixes the initial conditions
#
#     u(0)  = 1.071822516416917,
#     u'(0) = 1.710337353176786.
#
# This is the *equianharmonic* ℘ (lattice invariant g₂ = 0): a doubly
# periodic meromorphic function whose second-order poles sit on a
# rhombic lattice of equilateral triangles.  FW md:297 gives the real
# period: poles on the real axis at x = 1 + 2ωk, where
#
#     ω = Γ(1/3)³ / (2^(13/6) π) ≈ 1.363.
#
# The full equianharmonic lattice has primitive periods 2ω and
# 2ω·e^{iπ/3}, so the poles of u(z) = ℘(z + c₁) with c₁ = -1 sit at
#
#     z(m, n) = 1 + 2ω(m + n/2) + i · ω√3 · n,    m, n ∈ ℤ.
#
# (n = 0 recovers FW's real-axis row; n = ±1 the rows at Im z = ±ω√3
# ≈ ±2.36 visible in FW's figure.)  No pole-finding is needed — the
# lattice is analytic.
#
# The integrator path is `path_network_solve`'s visited tree for a
# single target at z = 30: the wedge walk threads the low-|u| passages
# *between* the pole rows, staying just below the real axis — exactly
# the path FW draws.  This is the same problem `test/pathnetwork_test.jl`
# pins quantitatively under PN.2.2 (u(30) to 2.13e-14 at BF-256).

using PadeTaylor
using CairoMakie
using Printf
using SpecialFunctions: gamma

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
const H_STEP = 0.5               # FW 2011 standard step length (md:279)
const ORDER  = 30                # FW 2011 standard Padé order (md:279)
const Z_END  = 30.0              # medium-range test target (md:299)
const U0     = 1.071822516416917 # ℘ IC, FW md:294
const UP0    = 1.710337353176786 # ℘ IC, FW md:295
const OUTPNG = joinpath(@__DIR__, "output", "fw2011_fig_5_1.png")

# Equianharmonic ℘ real half-period (FW md:297).
const Ω = gamma(1 / 3)^3 / (2.0^(13 / 6) * π)

# PI test problem with the z-term dropped: u'' = 6 u^2 (FW eq. 5.1).
pweier(z, u, up) = 6 * u^2

# ----------------------------------------------------------------------
# Analytic ℘ pole lattice over the displayed window
# ----------------------------------------------------------------------
# z(m, n) = 1 + 2ω(m + n/2) + i·ω√3·n.  Re z ∈ [-1, 31], Im z ∈ [-3, 3]
# ⇒ n ∈ {-1, 0, 1} (n = ±2 lands at Im ≈ ±4.7, off-window).
poles = ComplexF64[]
for n in -1:1, m in -2:18
    z = 1 + 2Ω * (m + n / 2) + im * (Ω * sqrt(3) * n)
    (-1.0 ≤ real(z) ≤ 31.0 && -3.0 ≤ imag(z) ≤ 3.0) && push!(poles, z)
end
@printf("FW 2011 Fig 5.1: %d ℘ poles in the window (real period 2ω = %.4f).\n",
        length(poles), 2Ω)

# ----------------------------------------------------------------------
# Integrator path: path_network_solve's visited tree to z = 30
# ----------------------------------------------------------------------
prob = PadeTaylorProblem(pweier, (U0, UP0), (0.0, Z_END); order = ORDER)
t0   = time()
# `path_network_solve` carries its own `order` kwarg (it does not read
# `prob.order`); pass it explicitly so the figure is order-30 by
# construction, not by coincidence with the kwarg default.
sol  = path_network_solve(prob, ComplexF64[Z_END + 0im];
                          h = H_STEP, order = ORDER,
                          max_steps_per_target = 400)
u30  = sol.grid_u[1]
relerr = abs(u30 - 1.095098255959744) / abs(1.095098255959744)
@printf("  path: %d visited nodes in %.2f s; u(30) rel-err %.2e (FW: 7.62e-14 at BF-256).\n",
        length(sol.visited_z), time() - t0, relerr)

# The visited tree for a single target is a single chain; plot it in
# visit order (visited_z[1] is the IC at the origin).
pathx = real.(sol.visited_z)
pathy = imag.(sol.visited_z)

# ----------------------------------------------------------------------
# Render
# ----------------------------------------------------------------------
fig = Figure(size = (1100, 320))
ax  = Axis(fig[1, 1]; xlabel = "Re z", ylabel = "Im z",
           title = "FW 2011 Fig. 5.1 — ℘ poles (dots) and the Padé-integrator path to z = 30",
           aspect = DataAspect(),
           limits = (-1, 31, -3, 3))
# Padé-integrator path: small dots joined in visit order.
lines!(ax, pathx, pathy; color = :black, linewidth = 0.9)
scatter!(ax, pathx, pathy; color = :black, markersize = 4)
# ℘ poles: large dots.
scatter!(ax, real.(poles), imag.(poles); color = :black, markersize = 15)

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
