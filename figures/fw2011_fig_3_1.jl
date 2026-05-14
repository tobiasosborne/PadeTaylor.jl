# figures/fw2011_fig_3_1.jl
#
# Reproduces Fornberg & Weideman 2011, Fig. 3.1 — the magnitude |u(z)|
# of the first Painlevé (PI) transcendent
#
#     u''(z) = 6 u(z)^2 + z
#
# for the near-tritronquée initial conditions
#
#     u(0)  = -0.1875,
#     u'(0) =  0.3049,
#
# displayed over the square z = x + i y, -10 <= x <= 10, -10 <= y <= 10.
#
# Source: references/markdown/FW2011_painleve_methodology_JCP230/
#         FW2011_painleve_methodology_JCP230.md:139-145.
#
# Because the ICs are CLOSE TO but not exactly the tritronquee values,
# FW 2011 reports pole fields in all five main sectors, with four of
# the five displaced some distance away from the origin (md:139).  The
# region immediately around the origin is comparatively pole-free.
# That qualitative signature — a calm central patch surrounded by five
# sectors of pole spikes, four of them set back from z = 0 — is the
# acceptance criterion (docs/figure_catalogue.md, FW 2011 Fig 3.1 row).
#
# Method.  We use `path_network_solve` (FW 2011 Section 3.1 two-stage
# path network) directly on the fine `N x N` lattice: Stage 1 builds
# the off-axis path tree that bridges the pole fields, Stage 2 takes a
# single Pade step from the nearest tree node up to each lattice
# point.  PI has real coefficients and the ICs sit on the real axis,
# so the solution obeys the Schwarz reflection u(zbar) = ubar(z); we
# pass `enforce_real_axis_symmetry = true`, which walks only the upper
# half-plane and mirrors — exact here, and roughly halves the work.
#
# `h = 0.5` and `order = 30` are FW 2011's working defaults (md:164,
# Section 5.1).  `N` is the lattice resolution; FW used 161.  It is the
# one knob worth turning — raise it for a finer surface at linear cost.

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
const U_CAP  = 50.0              # z-axis cap, matching FW's Fig 3.1 framing
const OUTPNG = joinpath(@__DIR__, "output", "fw2011_fig_3_1.png")

# PI right-hand side: u'' = 6 u^2 + z.
pI(z, u, up) = 6 * u^2 + z

# ----------------------------------------------------------------------
# Solve on the fine lattice
# ----------------------------------------------------------------------
xs, ys, grid = pole_field_lattice(HALF, N)

prob = PadeTaylorProblem(pI, (-0.1875, 0.3049), (0.0, 10.0); order = ORDER)

@printf("FW 2011 Fig 3.1: solving PI on a %d x %d lattice over [-%g, %g]^2 ...\n",
        N, N, HALF, HALF)
t0  = time()
sol = path_network_solve(prob, grid;
                         h = H_STEP,
                         max_steps_per_target = 4000,
                         enforce_real_axis_symmetry = true)
@printf("  done in %.1f s; %d visited tree nodes.\n",
        time() - t0, length(sol.visited_z))

# |u| surface, capped at U_CAP so the unbounded pole spikes do not
# flatten the rest of the field (this is exactly FW's Fig 3.1 framing).
absu = abs.(reshape(sol.grid_u, N, N))
ncov = count(isfinite, absu)
@printf("  Stage-2 coverage: %d / %d lattice points finite.\n", ncov, N * N)
Z = clamp.(absu, 0.0, U_CAP)

# Quantitative companion to the qualitative acceptance: the median
# |u| over a small central disc (|z| ≤ 2, the near-tritronquée
# pole-free patch) should be markedly below the median over an outer
# annulus (6 ≤ |z| ≤ 10, deep inside the displaced pole fields).
zc        = ComplexF64.(reshape(grid, N, N))
in_disc   = abs.(zc) .<= 2.0
in_annul  = 6.0 .<= abs.(zc) .<= 10.0
med_disc  = sort(vec(Z[in_disc]))[fld(count(in_disc), 2)]
med_annul = sort(vec(Z[in_annul]))[fld(count(in_annul), 2)]
@printf("  median |u|: central disc |z|≤2 = %.3f vs annulus 6≤|z|≤10 = %.3f (ratio %.1fx)\n",
        med_disc, med_annul, med_annul / med_disc)

# ----------------------------------------------------------------------
# Render (shared pole-field surface helper, figutil.jl)
# ----------------------------------------------------------------------
fig, _ = pole_field_figure(xs, ys, Z;
    title = "FW 2011 Fig. 3.1 — |u(z)| for PI, u(0) = -0.1875, u'(0) = 0.3049",
    ucap = U_CAP)

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
