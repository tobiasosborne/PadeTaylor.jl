# examples/tritronquee_3d.jl — 3D surface plot of the PI tritronquée
# pole field over a 2D complex window, in the style of FW 2011 Fig 3.1.
#
# Run with:
#   julia --project=examples examples/tritronquee_3d.jl
#
# Outputs:
#   examples/tritronquee_3d.png           — 3D surface log₁₀|u(z)|
#   examples/tritronquee_heatmap.png      — 2D heatmap log₁₀|u(z)|
#   examples/tritronquee_mask.png         — EdgeDetector pole-field mask
#
# Ground truth: PI tritronquée ICs per FW 2011 eq. (4.1)
# (`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:226`):
#   u(0)  = -0.1875543083404949
#   u'(0) =  0.3049055602612289
# ODE:  u'' = 6u² + z  (Painlevé I)
#
# This script uses the same setup as `test/phase9_tritronquee_test.jl`
# but on a denser grid (41×41 by default) for visual fidelity.

using PadeTaylor
using Plots

# --- Problem setup ---------------------------------------------------------

f_PI(z, u, up) = 6 * u^2 + z
u_tri  = -0.1875543083404949
up_tri =  0.3049055602612289

# Grid resolution.  EVEN N is intentional: with odd N the grid has a
# cell exactly at y = 0, and the walker computes that cell on a
# distinct path from cells at y = ±Δy.  In the pole-bearing region
# (positive real axis, where |u| varies wildly across adjacent
# cells), this manifests as a visible horizontal discontinuity at
# y = 0.  Even N avoids the problem entirely — the closest cells to
# the real axis sit at y = ±Δy/2, and the conjugate mirror handles
# them symmetrically with no special-case row.
# At [-20, 20]² (5× the linear extent of FW Fig 3.1, 25× the area),
# 500×500 gives FW Fig 3.1-style detail in ~40 s.  Drop to N = 250
# (~10 s) for faster iteration; N = 120 (~5 s) for a sanity check.
N = 500
xs = range(-20.0, 20.0; length = N)
ys = range(-20.0, 20.0; length = N)

# `zspan[2]` is unused for path-network-only callers; non-degenerate
# placeholder satisfies the PadeTaylorProblem guard.  Use the farthest
# grid corner.
zspan = (0.0 + 0.0im, ComplexF64(maximum(abs, xs) * sqrt(2)))
prob  = PadeTaylorProblem(f_PI, (u_tri, up_tri), zspan; order = 30)

# --- Solve with reflection-symmetry enforcement ---------------------------
# Real ICs (u₀, u'₀ ∈ ℝ) ⇒ u(z̄) = ū(z) globally (Schwarz reflection).
# The default path-network's `shuffle(rng, targets)` step creates an
# asymmetric walk-tree even when the problem is conjugate-symmetric
# (worklog 014 §"Bug 1"), and the asymmetric Stage-2 nearest-visited
# lookup then cascades 4-5 orders of magnitude into the per-cell |u|
# values for conjugate-pair cells.  The opt-in kwarg
# `enforce_real_axis_symmetry = true` (bead `padetaylor-dtj`) walks
# only upper-half + on-axis targets and mirrors lower-half via conj —
# bit-exact symmetric, roughly 2× faster than a full-plane walk, and
# correct *because* the ODE coefficients are real.  For non-real ICs /
# complex-coefficient ODEs the kwarg is inapplicable; the call would
# throw on the IC validation, so leave it `false`.
#
# Even N is still recommended.  With odd N the grid has a cell exactly
# at y = 0; the wedge walker computes that cell on a special direction
# (the goal-aligned axis), producing the "low-level ridges in flat
# areas" artifact (FW2011...md:208, worklog 014 §"Bug 2") visible as a
# horizontal stripe of false-pole flags in the asymptotic smooth region.
# Even N avoids the y = 0 row entirely.

@assert iseven(N) "N must be even — see the comment above (worklog 014 §Bug 2)."
grid_mat = ComplexF64[xs[i] + im * ys[j] for i in 1:N, j in 1:N]
grid_vec = vec(grid_mat)   # column-major; reshape preserves the matrix layout

println("Solving PI tritronquée on $(N)×$(N) grid over [$(Int(minimum(xs))), $(Int(maximum(xs)))]² (enforce_real_axis_symmetry)...")
@time sol = path_network_solve(prob, grid_vec; h = 0.5,
                                max_steps_per_target = 2000,
                                enforce_real_axis_symmetry = true)
println("  Total targets: $(length(grid_vec))")
println("  Visited nodes (upper-half only): $(length(sol.visited_z))")
println("  Grid coverage: $(count(isfinite.(real.(sol.grid_u))))/$(length(grid_vec))")

u_grid = reshape(sol.grid_u, (N, N))
h_grid = step(xs)

# Clamp the pole spikes so the surface plot doesn't get dominated by a
# single mountain.  log₁₀|u(z)| is the natural display — FW Fig 3.3
# convention.  Clip at +6 to limit the dynamic range.
log_abs_u = clamp.(log10.(abs.(u_grid)), -1.0, 6.0)

# Pole-field mask via EdgeDetector (FW eq. 3.3 + Fig 3.3 contour at
# level 0.001 on log₁₀|Δu|).
mask = pole_field_mask(u_grid, h_grid)

# --- 3D surface plot -------------------------------------------------------

println("Rendering 3D surface (log₁₀|u(z)|)...")
# Plots stores matrices as z[y, x] for surface(x, y, z), so transpose.
plt_3d = surface(
    xs, ys, transpose(log_abs_u);
    xlabel = "x = Re(z)",
    ylabel = "y = Im(z)",
    zlabel = "log₁₀|u(z)|",
    title  = "PI tritronquée — u'' = 6u² + z, IC (FW eq. 4.1)\nlog₁₀|u(z)| over [$(Int(minimum(xs))), $(Int(maximum(xs)))]²",
    color  = :viridis,
    camera = (40, 35),       # azimuth, elevation
    size   = (1000, 800),
    clims  = (-1, 6),
)
savefig(plt_3d, "examples/tritronquee_3d.png")
println("  Wrote examples/tritronquee_3d.png")

# --- 2D heatmap -------------------------------------------------------------

println("Rendering 2D heatmap...")
plt_heat = heatmap(
    xs, ys, transpose(log_abs_u);
    xlabel = "x = Re(z)",
    ylabel = "y = Im(z)",
    title  = "PI tritronquée — log₁₀|u(z)|",
    color  = :viridis,
    aspect_ratio = :equal,
    size   = (800, 700),
    clims  = (-1, 6),
)
savefig(plt_heat, "examples/tritronquee_heatmap.png")
println("  Wrote examples/tritronquee_heatmap.png")

# --- Pole-field mask --------------------------------------------------------

println("Rendering EdgeDetector pole-field mask...")
plt_mask = heatmap(
    xs, ys, transpose(Float64.(mask));
    xlabel = "x = Re(z)",
    ylabel = "y = Im(z)",
    title  = "PI tritronquée — EdgeDetector mask (FW eq. 3.3, level=0.001)",
    color  = :greys,
    aspect_ratio = :equal,
    size   = (800, 700),
    legend = false,
)
savefig(plt_mask, "examples/tritronquee_mask.png")
println("  Wrote examples/tritronquee_mask.png")

println("\nDone. View the PNGs from your host filesystem.")
