# examples/tritronquee_bvp_compose.jl — same PI tritronquée problem
# as `tritronquee_3d.jl`, but with per-row BVP fill applied to interior
# smooth runs after the upper-half path-network walk.  Demonstrates FW
# 2011 §3.2 + §4.4's prescribed cure for the "low-level ridges in flat
# areas" artifact (FW2011...md:208).
#
# Renders three side-by-side comparisons:
#   examples/tritronquee_before_bvp.png — bare IVP (path-network only)
#   examples/tritronquee_after_bvp.png  — IVP + interior BVP fill
#   examples/tritronquee_bvp_mask.png   — EdgeDetector mask (unchanged)
#
# Wall ~30 s at N = 250 over [-20, 20]².

using PadeTaylor
using PadeTaylor.PathNetwork: path_network_solve
using PadeTaylor.EdgeDetector: pole_field_mask
using PadeTaylor.BVP: bvp_solve
using Plots

# --- Problem -------------------------------------------------------------

f_PI(z, u, up) = 6 * u^2 + z       # 2nd-order form for IVP path-network
f_bvp(z, u)    = 6 * u^2 + z       # 1st-order form for BVP (autonomous in u')
∂f_bvp(z, u)   = 12 * u

u_tri  = -0.1875543083404949
up_tri =  0.3049055602612289

# Use even N to avoid the y=0 row issue documented in tritronquee_3d.jl.
N  = 250
xs = range(-20.0, 20.0; length = N)
ys = range(-20.0, 20.0; length = N)
zspan = (0.0 + 0.0im, ComplexF64(20 * sqrt(2)))
prob  = PadeTaylorProblem(f_PI, (u_tri, up_tri), zspan; order = 30)

# --- IVP path-network (upper half + mirror) ------------------------------

upper_grid  = ComplexF64[]
upper_to_ij = Tuple{Int, Int}[]
for i in 1:N, j in 1:N
    if ys[j] > 0
        push!(upper_grid, xs[i] + im * ys[j])
        push!(upper_to_ij, (i, j))
    end
end

println("IVP solve: $(N)×$(N) over [-20, 20]², upper-half walk + mirror...")
@time sol = path_network_solve(prob, upper_grid; h = 0.5,
                                max_steps_per_target = 2000)
println("  Visited: $(length(sol.visited_z)), grid coverage: $(count(isfinite.(real.(sol.grid_u))))/$(length(upper_grid))")

u_grid = Matrix{ComplexF64}(undef, N, N)
for (idx, (i, j)) in enumerate(upper_to_ij)
    u_grid[i, j]         = sol.grid_u[idx]
    u_grid[i, N + 1 - j] = conj(sol.grid_u[idx])
end
h_grid = step(xs)

# Snapshot the bare-IVP field for the "before" plot.
u_grid_ivp = copy(u_grid)
mask = pole_field_mask(u_grid, h_grid)
println("Mask cells (IVP, before BVP fill): $(count(mask)) / $(N^2)")

# --- Per-row BVP fill on interior smooth runs ---------------------------

# A "bridgeable" smooth run in row j is a maximal contiguous block of
# mask=false cells flanked by mask=true cells on BOTH sides, with both
# flanks strictly interior to the grid (so the EdgeDetector mask there
# is well-defined).  For each such run, run bvp_solve with Dirichlet
# BCs from the IVP values at the flanks, and replace the interior
# u-values with the BVP's barycentric interpolation at each cell.
println("Applying per-row BVP fill...")
bvp_count = 0
@time for j in 2:(N - 1)
    i = 2
    while i ≤ N - 1
        if mask[i, j]
            i += 1
            continue
        end
        run_start = i
        while i ≤ N - 1 && !mask[i, j]
            i += 1
        end
        run_end = i - 1
        left_flank  = run_start - 1
        right_flank = run_end + 1
        if 2 ≤ left_flank ≤ N - 1 && 2 ≤ right_flank ≤ N - 1 &&
           mask[left_flank, j] && mask[right_flank, j]

            z_a = ComplexF64(xs[left_flank]  + im * ys[j])
            z_b = ComplexF64(xs[right_flank] + im * ys[j])
            u_a = u_grid[left_flank,  j]
            u_b = u_grid[right_flank, j]
            try
                bvp_sol = bvp_solve(f_bvp, ∂f_bvp, z_a, z_b, u_a, u_b; N = 20)
                for k in run_start:run_end
                    z_k = ComplexF64(xs[k] + im * ys[j])
                    u_k, _ = bvp_sol(z_k)
                    u_grid[k, j] = u_k
                end
                bvp_count += 1
            catch err
                # A few BVPs may fail to converge near the pole-edge cells
                # where the IVP-derived BC values are themselves near-pole
                # magnitudes; tag them as "skipped" and leave the IVP value.
                # In production lattice_dispatch_solve such runs would be
                # flagged :ivp_only.
            end
        end
    end
end
println("  BVP solves: $bvp_count")

mask_after = pole_field_mask(u_grid, h_grid)
println("Mask cells (after BVP fill): $(count(mask_after)) / $(N^2)")

# --- Render --------------------------------------------------------------

function render(u, fname; title_extra = "")
    log_abs_u = clamp.(log10.(abs.(u)), -1.0, 6.0)
    plt = heatmap(xs, ys, transpose(log_abs_u);
        xlabel = "x = Re(z)",
        ylabel = "y = Im(z)",
        title  = "PI tritronquée — log₁₀|u(z)|$title_extra",
        color  = :viridis,
        aspect_ratio = :equal,
        size   = (800, 700),
        clims  = (-1, 6),
    )
    savefig(plt, fname)
    println("  Wrote $fname")
end

println("Rendering...")
render(u_grid_ivp, "examples/tritronquee_before_bvp.png";
       title_extra = " (path-network only)")
render(u_grid, "examples/tritronquee_after_bvp.png";
       title_extra = " (path-network + interior BVP fill)")

# Also render the mask after BVP fill for direct comparison
plt_mask = heatmap(xs, ys, transpose(Float64.(mask_after));
    xlabel = "x = Re(z)", ylabel = "y = Im(z)",
    title  = "EdgeDetector mask (after BVP fill)",
    color  = :greys, aspect_ratio = :equal, size = (800, 700), legend = false)
savefig(plt_mask, "examples/tritronquee_bvp_mask.png")
println("  Wrote examples/tritronquee_bvp_mask.png")

println("Done.")
