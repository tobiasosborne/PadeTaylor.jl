# external/probes/loop-closure-fig1/smoke.jl
#
# Tiny sanity check for the three architectural unknowns flagged in the
# probe plan.  Designed to crash loudly on any unmet assumption.

using PadeTaylor
using PadeTaylor.CoordTransforms: pIII_transformed_rhs, pIII_z_to_ζ
using DelaunayTriangulation
using Printf

const α, β, γ, δ = -0.5, -0.5, 1.0, -1.0
const Z0, U0, UP0 = 1.0 + 0.0im, 0.25 + 0.0im, 1.0 + 0.0im
const ORDER = 30

# Minimal solve: tiny target list — we just need a small visited tree to
# poke the data structures.
ζ0, w0, wp0 = pIII_z_to_ζ(Z0, U0, UP0)
@printf("IC ζ-frame: (ζ, w, w') = (%s, %s, %s)\n", ζ0, w0, wp0)

f = pIII_transformed_rhs(α, β, γ, δ)
prob = PadeTaylorProblem(f, (w0, wp0), (ζ0, complex(2.0, 1.0)); order = ORDER)
R(ζ) = max(0.1, (8.0 - real(ζ)) / 20.0)

# Walk to a few targets, just enough to populate the visited tree.
small_targets = ComplexF64[complex(x, y) for y in -1.0:0.5:1.0
                                          for x in -1.0:0.5:1.5]
sol = path_network_solve(prob, small_targets;
                          h = R(ζ0),
                          node_separation = R,
                          step_size_policy = :adaptive_ffw,
                          adaptive_tol = 1e-10,
                          k_conservative = 1e-3,
                          max_rescales = 50,
                          max_steps_per_target = 8000)

@printf("Visited tree: %d nodes\n", length(sol.visited_z))
@printf("sol.visited_sheet isa %s\n", typeof(sol.visited_sheet))
@printf("length(sol.visited_sheet) = %d\n", length(sol.visited_sheet))
@printf("sol.visited_sheet[1] = %s\n", sol.visited_sheet[1])
@printf("length(sol.visited_sheet[1]) == 0 ? %s\n",
        length(sol.visited_sheet[1]) == 0)

# Unknown 2: can we call _evaluate_pade directly?
P = sol.visited_pade[1]
hv = sol.visited_h[1]
t = (complex(0.05, 0.0) - sol.visited_z[1]) / hv
u_test = PadeTaylor.PathNetwork._evaluate_pade(P, t)
@printf("_evaluate_pade(visited_pade[1], t) = %s\n", u_test)

# Unknown 3: DelaunayTriangulation.triangulate API smoke.
pts = [(real(z), imag(z)) for z in sol.visited_z]
tri = triangulate(pts)
@printf("DelaunayTriangulation.triangulate built; type = %s\n",
        typeof(tri))
edges_set = each_edge(tri)
n_edges = length(collect(edges_set))
@printf("Unique edges from each_edge(): %d (visited n=%d ⇒ expect ~3n)\n",
        n_edges, length(sol.visited_z))

# Pull one edge and check its tuple shape.
some_edge = first(edges_set)
@printf("Sample edge: %s (type %s)\n", some_edge, typeof(some_edge))
println("SMOKE OK")
