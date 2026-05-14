# figures/fw2011_fig_4_3.jl
#
# Reproduces Fornberg & Weideman 2011, Fig. 4.3 — the NIST Handbook
# example: the magnitude |u(z)| of the PI transcendent
#
#     u''(z) = 6 u(z)^2 + z
#
# for the near-tronquée initial conditions
#
#     u(0)  = 0,
#     u'(0) = 1.8518,
#
# over the square z = x + i y, -10 <= x,y <= 10.
#
# Source: references/markdown/FW2011_painleve_methodology_JCP230/
#         FW2011_painleve_methodology_JCP230.md:231-245.
#
# This is one of the two NIST Handbook tronquée-family examples (FW
# md:231-233, "similar to Figs. 32.3.3 and 32.3.4 in [21]").  The ICs
# u'(0) = 1.8518 sit just below the exact tronquée value
# u'(0) = 1.85185403375822 — see `fw2011_fig_4_2.jl`.  As u'(0)
# approaches the tronquée case from below, the left pole field is
# pushed far out along the negative real axis and returns with a row
# of poles ON the negative real axis (FW md:233); the right pole
# field is essentially unchanged.  Fig. 4.4 is the same display for
# u'(0) = 1.8519 — just ABOVE the tronquée value — and FW notes the
# pole count there differs from this figure by exactly one.
#
# The two annotations FW overlays (md:233): a thick line along the
# real axis (the |u(x)| profile — this is the dashed curve of
# Fig. 4.2a traced on the surface) and a short cross-line on the
# imaginary axis marking the origin.
#
# Method: `path_network_solve` on the fine lattice with Schwarz
# symmetry (PI has real coefficients, the ICs are on the real axis),
# `h = 0.5`, `order = 30` — FW 2011's working defaults.

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
const U_CAP  = 50.0              # z-axis cap, matching FW's Fig 4.3 framing
const UP0    = 1.8518            # near-tronquée IC, FW Fig 4.3 (md:245)
const OUTPNG = joinpath(@__DIR__, "output", "fw2011_fig_4_3.png")

# PI right-hand side: u'' = 6 u^2 + z.
pI(z, u, up) = 6 * u^2 + z

# ----------------------------------------------------------------------
# Solve on the fine lattice
# ----------------------------------------------------------------------
xs, ys, grid = pole_field_lattice(HALF, N)
prob = PadeTaylorProblem(pI, (0.0, UP0), (0.0, 10.0); order = ORDER)

@printf("FW 2011 Fig 4.3: solving PI (u(0)=0, u'(0)=%g) on a %d x %d lattice ...\n",
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
    title = "FW 2011 Fig. 4.3 — |u(z)| for PI, u(0) = 0, u'(0) = $(UP0)",
    ucap = U_CAP)

# Thick real-axis line: |u(x)| traced on the surface (y = 0 is lattice
# row j with ys[j] = 0; N is odd so this is the exact midpoint).
jzero = (N + 1) ÷ 2
@assert ys[jzero] == 0.0
lines!(ax, collect(xs), zeros(N), Z[:, jzero]; color = :black, linewidth = 3)
# Short cross-line on the imaginary axis marking the origin.
lines!(ax, [0.0, 0.0], [-1.0, 1.0], [0.0, 0.0]; color = :black, linewidth = 3)

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
