# figures/fw2011_fig_4_4.jl
#
# Reproduces Fornberg & Weideman 2011, Fig. 4.4 — the companion to
# Fig. 4.3.  Same display, same equation
#
#     u''(z) = 6 u(z)^2 + z,
#
# same NIST Handbook setting (u(0) = 0), but the initial slope is now
#
#     u'(0) = 1.8519
#
# — just ABOVE the exact tronquée value u'(0) = 1.85185403375822,
# whereas Fig. 4.3's u'(0) = 1.8518 is just below it.
#
# Source: references/markdown/FW2011_painleve_methodology_JCP230/
#         FW2011_painleve_methodology_JCP230.md:253-255 (+ :231-237).
#
# Why the pair matters.  The tronquée case is a sharp transition: as
# u'(0) crosses 1.85185403375822, the left pole field's front row of
# poles flips between "a gap straddling the negative real axis" and
# "a pole sitting on the negative real axis" (FW md:233, md:259 Fig.
# 4.5 caption).  The catalogue acceptance for this figure
# (`docs/figure_catalogue.md`, FW 2011 Fig 4.4 row) is precisely that
# the pole count over the displayed window differs from Fig. 4.3 by
# *exactly one* — the single extra/missing front-row pole.
#
# The annotations and method are identical to `fw2011_fig_4_3.jl`;
# only `UP0` changes.  Run both and compare the printed pole proxies
# and the rendered surfaces side by side.

using PadeTaylor
using CairoMakie
using Printf
include(joinpath(@__DIR__, "figutil.jl"))

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
const N      = 121               # lattice is N x N over [-10, 10]^2 (FW: 161)
const HALF   = 10.0
const H_STEP = 0.5               # FW 2011 working step length (md:164)
const ORDER  = 30                # FW 2011 Section 5.1 Taylor order
const U_CAP  = 50.0              # z-axis cap, matching FW's Fig 4.4 framing
const UP0    = 1.8519            # near-tronquée IC, FW Fig 4.4 (md:255)
const OUTPNG = joinpath(@__DIR__, "output", "fw2011_fig_4_4.png")

# PI right-hand side: u'' = 6 u^2 + z.
pI(z, u, up) = 6 * u^2 + z

# ----------------------------------------------------------------------
# Solve on the fine lattice
# ----------------------------------------------------------------------
xs, ys, grid = pole_field_lattice(HALF, N)
prob = PadeTaylorProblem(pI, (0.0, UP0), (0.0, 10.0); order = ORDER)

@printf("FW 2011 Fig 4.4: solving PI (u(0)=0, u'(0)=%g) on a %d x %d lattice ...\n",
        UP0, N, N)
t0  = time()
sol = path_network_solve(prob, grid;
                         h = H_STEP,
                         max_steps_per_target = 4000,
                         enforce_real_axis_symmetry = true)
@printf("  done in %.1f s; %d visited tree nodes.\n",
        time() - t0, length(sol.visited_z))

absu = abs.(reshape(sol.grid_u, N, N))           # absu[i, j] at (xs[i], ys[j])
ncov = count(isfinite, absu)
ncap = count(>=(U_CAP), filter(isfinite, vec(absu)))
@printf("  Stage-2 coverage: %d / %d finite; %d cells at/above the |u|=%g cap (pole proxy).\n",
        ncov, N * N, ncap, U_CAP)
Z = clamp.(absu, 0.0, U_CAP)

# ----------------------------------------------------------------------
# Render (shared pole-field surface helper, figutil.jl) + FW's two
# real-/imaginary-axis annotations.
# ----------------------------------------------------------------------
fig, ax = pole_field_figure(xs, ys, Z;
    title = "FW 2011 Fig. 4.4 — |u(z)| for PI, u(0) = 0, u'(0) = $(UP0)",
    ucap = U_CAP)

# Thick real-axis line: |u(x)| traced on the surface.
jzero = (N + 1) ÷ 2
@assert ys[jzero] == 0.0
lines!(ax, collect(xs), zeros(N), Z[:, jzero]; color = :black, linewidth = 3)
# Short cross-line on the imaginary axis marking the origin.
lines!(ax, [0.0, 0.0], [-1.0, 1.0], [0.0, 0.0]; color = :black, linewidth = 3)

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
