# figures/fw2011_fig_3_2.jl
#
# Reproduces Fornberg & Weideman 2011, Fig. 3.2 — the Stage-1
# integration-path tree of the FW 2011 Section 3.1 path network.
#
# Source: references/markdown/FW2011_painleve_methodology_JCP230/
#         FW2011_painleve_methodology_JCP230.md:151-164.
#
# The same PI problem as Fig. 3.1 (u'' = 6 u^2 + z, near-tritronquee
# ICs u(0) = -0.1875, u'(0) = 0.3049) is used, but here the interest
# is the path TOPOLOGY rather than the solution magnitude.  FW 2011
# Stage 1 (md:155-164): lay a coarse 40 x 40 grid across [-10, 10]^2,
# and from the origin grow a tree of length-`h` steps that reaches
# within `h` of every coarse node.  Each step picks, among five wedge
# directions (straight at the goal, +/-22.5 deg, +/-45 deg), the one
# minimising |u| — the heuristic for "stay in the low-magnitude
# passages between poles".  A new target's walk starts from the
# nearest already-visited node, so the paths form a single tree rooted
# at the origin and never cross (md:164).
#
# FW 2011 uses h = 0.3 for this illustration specifically (md:158,
# md:164) and reports "about 1540 steps were needed for 1600 target
# points".  That step count is the quantitative sanity check; the
# acceptance proper is the visual topology match — a connected,
# non-crossing tree rooted at the origin (docs/figure_catalogue.md,
# FW 2011 Fig 3.2 row).
#
# The tree edges come from `PathNetworkSolution.visited_parent`, added
# to the solver specifically so this figure can be drawn: edge k is
# the segment from `visited_z[visited_parent[k]]` to `visited_z[k]`.

using PadeTaylor
using CairoMakie
using Printf

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
const N      = 40                # FW 2011's coarse grid is exactly 40 x 40
const HALF   = 10.0
const H_STEP = 0.3               # FW 2011 Fig 3.2 step length (md:158, md:164)
const ORDER  = 30                # FW 2011 Section 5.1 Taylor order
const OUTPNG = joinpath(@__DIR__, "output", "fw2011_fig_3_2.png")

# PI right-hand side: u'' = 6 u^2 + z.
pI(z, u, up) = 6 * u^2 + z

# ----------------------------------------------------------------------
# Build the Stage-1 path tree over the coarse grid
# ----------------------------------------------------------------------
xs = range(-HALF, HALF; length = N)
ys = range(-HALF, HALF; length = N)
grid = ComplexF64[complex(x, y) for y in ys for x in xs]

prob = PadeTaylorProblem(pI, (-0.1875, 0.3049), (0.0, 10.0); order = ORDER)

@printf("FW 2011 Fig 3.2: growing the Stage-1 path tree over a %d x %d grid ...\n",
        N, N)
t0  = time()
# Plain FW 2011 algorithm: full-plane random-target walk (NOT the
# Schwarz-symmetry shortcut) — Fig 3.2 illustrates the actual
# asymmetric tree the paper's algorithm produces.
sol = path_network_solve(prob, grid;
                         h = H_STEP,
                         max_steps_per_target = 4000)
n_nodes = length(sol.visited_z)
n_steps = n_nodes - 1            # every node bar the root is one step
@printf("  done in %.1f s; %d steps for %d target points (FW 2011: ~1540 for 1600).\n",
        time() - t0, n_steps, N * N)

# ----------------------------------------------------------------------
# Render: coarse grid as dots, the path tree as line segments, the
# origin (root) as a solid black circle.
# ----------------------------------------------------------------------
# Tree edges: for each non-root node k, the segment parent -> k.
seg = Point2f[]
for k in 2:n_nodes
    p = sol.visited_parent[k]
    push!(seg, Point2f(real(sol.visited_z[p]), imag(sol.visited_z[p])))
    push!(seg, Point2f(real(sol.visited_z[k]), imag(sol.visited_z[k])))
end

fig = Figure(size = (760, 760))
ax  = Axis(fig[1, 1];
           xlabel = "Re z", ylabel = "Im z", aspect = DataAspect(),
           title  = "FW 2011 Fig. 3.2 — Stage-1 path tree, h = 0.3, $(N)x$(N) grid",
           limits = (-HALF - 0.4, HALF + 0.4, -HALF - 0.4, HALF + 0.4))
scatter!(ax, real.(grid), imag.(grid); color = :black, markersize = 2.5)
linesegments!(ax, seg; color = :black, linewidth = 0.6)
scatter!(ax, [0.0], [0.0]; color = :black, markersize = 11)

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
