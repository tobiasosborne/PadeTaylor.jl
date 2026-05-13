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

# Grid resolution.  At [-20, 20]² (5× the linear extent of FW Fig 3.1,
# 25× the area), 501×501 (h_grid = 0.08) gives FW Fig 3.1-style detail
# (their figure was 161×161 over [-10, 10]²; we're a hair denser per
# unit area).  Wall ≈ 100 s.  Drop to N = 251 (~25 s) for faster
# iteration; N = 121 (~11 s) for a quick sanity check.
N = 501
xs = range(-20.0, 20.0; length = N)
ys = range(-20.0, 20.0; length = N)

# `zspan[2]` is unused for path-network-only callers; non-degenerate
# placeholder satisfies the PadeTaylorProblem guard.  Use the farthest
# grid corner.
zspan = (0.0 + 0.0im, ComplexF64(maximum(abs, xs) * sqrt(2)))
prob  = PadeTaylorProblem(f_PI, (u_tri, up_tri), zspan; order = 30)

# --- Solve with reflection-symmetry enforcement ---------------------------
# Real ICs (u₀, u'₀ ∈ ℝ) ⇒ u(z̄) = ū(z) globally (Schwarz reflection).
# The path-network's `shuffle(rng, targets)` step creates an asymmetric
# walk-tree even when the problem is conjugate-symmetric (e.g., the
# tree gains ~25 more nodes on one half-plane than the other for a
# 81×81 [-20, 20]² PI tritronquée grid), and the asymmetric Stage-2
# nearest-visited lookup then cascades 4-5 orders of magnitude into
# the per-cell |u| values for conjugate-pair cells.  Fix: walk only the
# upper half + real axis; mirror to the lower half via `u(z̄) = ū(z)`.
# This is bit-exact symmetric, runs roughly 2× faster than walking the
# full grid, and is correct *because* the ODE coefficients are real.
# (For non-real ICs / complex-coefficient ODEs the mirror trick is
# inapplicable; use `path_network_solve` directly on the full grid.)

upper_grid  = ComplexF64[]
upper_to_ij = Tuple{Int, Int}[]
for i in 1:N, j in 1:N
    if ys[j] > 0   # strictly upper half — y=0 row reconstructed from
                   # Re(u(x, +Δy)) below, NOT walked independently.
        push!(upper_grid, xs[i] + im * ys[j])
        push!(upper_to_ij, (i, j))
    end
end

println("Solving PI tritronquée on $(N)×$(N) grid over [$(Int(minimum(xs))), $(Int(maximum(xs)))]² (upper-half walk + conjugate mirror + y=0 reconstruction)...")
@time sol = path_network_solve(prob, upper_grid; h = 0.5,
                                max_steps_per_target = 2000)
println("  Upper-half targets: $(length(upper_grid))")
println("  Visited nodes: $(length(sol.visited_z))")
println("  Grid coverage: $(count(isfinite.(real.(sol.grid_u))))/$(length(upper_grid))")

u_grid = Matrix{ComplexF64}(undef, N, N)
for (idx, (i, j)) in enumerate(upper_to_ij)
    u_grid[i, j] = sol.grid_u[idx]
    u_grid[i, N + 1 - j] = conj(sol.grid_u[idx])   # mirror to lower
end

# Reconstruct the y=0 row by Schwarz reflection.  Real ICs ⇒ u(x) ∈ ℝ
# for x ∈ ℝ; Taylor expansion off the real axis: u(x + iΔy) = u(x) +
# i·Δy·u'(x) + O(Δy²).  So u(x, 0) ≈ Re(u(x, +Δy)).  Walking y=0
# directly produces a path-dependent discontinuity vs y = ±Δy
# (different visited-tree branches near the singular pole-free
# direction); reconstructing from the +Δy row keeps the y=0 line
# continuous with its neighbours at O(Δy²) cost.
j_axis  = findfirst(==(0.0), ys)
if j_axis !== nothing
    for i in 1:N
        u_grid[i, j_axis] = complex(real(u_grid[i, j_axis + 1]), 0.0)
    end
end
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
