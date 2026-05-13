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

# Grid resolution.  41×41 takes ~5 s on this hardware; 51×51 ~10 s.
N = 41
xs = range(-4.0, 4.0; length = N)
ys = range(-4.0, 4.0; length = N)

# Flat-vector form for the path-network driver.
grid_mat = ComplexF64[xs[i] + im * ys[j] for i in 1:N, j in 1:N]
grid_vec = vec(grid_mat)

# `zspan[2]` is unused for path-network-only callers; non-degenerate
# placeholder satisfies the PadeTaylorProblem guard.
zspan = (0.0 + 0.0im, ComplexF64(4 * sqrt(2)))
prob  = PadeTaylorProblem(f_PI, (u_tri, up_tri), zspan; order = 30)

# --- Solve -----------------------------------------------------------------

println("Solving PI tritronquée on $(N)×$(N) grid over [-4, 4]²...")
@time sol = path_network_solve(prob, grid_vec; h = 0.5)
println("  Visited nodes: $(length(sol.visited_z))")
println("  Grid coverage: $(count(isfinite.(real.(sol.grid_u))))/$(N*N)")

# Reshape to 2D matrix; u_grid[i, j] = u(xs[i] + im·ys[j]).
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
    title  = "PI tritronquée — u'' = 6u² + z, IC (FW eq. 4.1)\nlog₁₀|u(z)| over [-4, 4]²",
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
