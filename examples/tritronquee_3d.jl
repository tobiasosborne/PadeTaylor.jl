# examples/tritronquee_3d.jl — the PI tritronquée pole field over a 2D
# complex window, in the style of FW 2011 Fig 3.1, computed with the
# edge-gated pole-field solver.
#
# Run with:
#   julia --project=examples examples/tritronquee_3d.jl
#
# Outputs:
#   examples/tritronquee_3d.png           — 3D surface log₁₀|u(z)|
#   examples/tritronquee_heatmap.png      — 2D heatmap log₁₀|u(z)|
#   examples/tritronquee_mask.png         — the region-grown pole-field mask
#
# Ground truth: PI tritronquée ICs per FW 2011 eq. (4.1)
# (`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:226`):
#   u(0)  = -0.1875543083404949
#   u'(0) =  0.3049055602612289
# ODE:  u'' = 6u² + z  (Painlevé I)
#
# ----------------------------------------------------------------------
# Why this script uses `edge_gated_pole_field_solve`, not the plain walk
#
# An earlier version of this script solved the *whole* window with a
# bare `path_network_solve`.  That left a visible artifact: faint
# thread-like seams of slightly-elevated |u| running through the
# pole-free sectors — the "low-level ridges in flat areas" that FW 2011
# warns about (md:208), and that the EdgeDetector mask flagged as
# spurious one-cell-wide false pole fields.  They are walk-tree seams:
# where two regions of the Stage-2 nearest-visited lookup meet, the
# per-cell |u| jumps, because the IVP path-network is *unstable in
# smooth regions* (md:401) and accumulates error there.
#
# `EdgeGatedSolve.edge_gated_pole_field_solve` is FW's own cure: it
# grows the IVP's target set outward from the initial condition, one
# frontier at a time, admitting only cells the §3.2.2 edge detector
# confirms as pole-field — and only cells connected to the field
# already found.  The IVP therefore never targets a smooth cell, the
# seams never form, and `egs.field_mask` is the clean, region-grown,
# flood-filled pole field (no morphological-opening specks, no threads).
#
# The trade-off is deliberate and documented in the EdgeGatedSolve
# module docstring: the pole-free sectors are *not solved at all*
# (`egs.u_grid` is `NaN` there).  For this picture that is correct —
# the tritronquée is pole-free in four of its five sectors, and a
# headline |u| surface only needs to show *where the poles are*.  The
# unsolved sectors are rendered at the display floor (deep purple),
# the same display convention the pole-spike clamp at +6 already uses.

using PadeTaylor
using Plots

# --- Problem setup ---------------------------------------------------------

f_PI(z, u, up) = 6 * u^2 + z
u_tri  = -0.1875543083404949
up_tri =  0.3049055602612289

# Grid resolution.  The §3.2.2 edge detector needs a lattice spacing
# of roughly ≲ 1 to separate pole-field from smooth (see the
# EdgeGatedSolve "Grid resolution matters" docstring); finer gives
# more detail but the region-growing solve costs more per pass.  At
# [-20, 20]² (5× the linear extent of FW Fig 3.1), N = 161 gives a
# spacing of 0.25 and renders in a couple of minutes.  Drop to N = 81
# (spacing 0.5) for a faster sanity pass.
N = 161
xs = range(-20.0, 20.0; length = N)
ys = range(-20.0, 20.0; length = N)

# `zspan[2]` is unused for path-network-only callers; the IC at
# `zspan[1]` is what the edge-gated seed grows out from.
zspan = (0.0 + 0.0im, ComplexF64(maximum(abs, xs) * sqrt(2)))
prob  = PadeTaylorProblem(f_PI, (u_tri, up_tri), zspan; order = 30)

# --- Edge-gated pole-field solve ------------------------------------------
# Region-growing from a seed at the IC; the IVP path-network is confined
# to edge-detector-confirmed pole-field cells, so the smooth-region
# walk-tree seams never form.  `grow_rings = 4` is the FW-Fig-4.7
# working width.

println("Edge-gated pole-field solve: PI tritronquée on $(N)×$(N) grid " *
        "over [$(Int(minimum(xs))), $(Int(maximum(xs)))]²...")
@time egs = edge_gated_pole_field_solve(prob, xs, ys;
                                        h = 0.5, order = 30,
                                        grow_rings = 4, verbose = true)
println("  Region-growing passes: $(egs.iterations)")
println("  Pole-field cells:      $(count(egs.field_mask)) / $(N^2)")
println("  Visited nodes:         $(length(egs.pn_solution.visited_z))")

u_grid = egs.u_grid       # NaN in the unsolved (pole-free) sectors
mask   = egs.field_mask   # the clean, region-grown pole-field mask

# log₁₀|u(z)| is the natural display (FW Fig 3.3 convention).  Clamp
# the pole spikes at +6 so a single mountain does not dominate; render
# the unsolved pole-free sectors at the display floor −1 (see the
# header comment — those sectors are pole-free by construction).
log_abs_u = map(u_grid) do u
    isnan(u) ? -1.0 : clamp(log10(abs(u)), -1.0, 6.0)
end

# --- 3D surface plot -------------------------------------------------------

println("Rendering 3D surface (log₁₀|u(z)|)...")
# Plots stores matrices as z[y, x] for surface(x, y, z), so transpose.
plt_3d = surface(
    xs, ys, transpose(log_abs_u);
    xlabel = "x = Re(z)",
    ylabel = "y = Im(z)",
    zlabel = "log₁₀|u(z)|",
    title  = "PI tritronquée — u'' = 6u² + z, IC (FW eq. 4.1)\n" *
             "log₁₀|u(z)| over [$(Int(minimum(xs))), $(Int(maximum(xs)))]² " *
             "(edge-gated pole-field solve)",
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
# `egs.field_mask` is the region-grown, flood-filled pole field — the
# clean counterpart of the raw `pole_field_mask` edge classifier, with
# the spurious smooth-region threads already removed by the
# morphological opening + flood-fill inside `edge_gated_pole_field_solve`.

println("Rendering region-grown pole-field mask...")
plt_mask = heatmap(
    xs, ys, transpose(Float64.(mask));
    xlabel = "x = Re(z)",
    ylabel = "y = Im(z)",
    title  = "PI tritronquée — region-grown pole-field mask (edge-gated solve)",
    color  = :greys,
    aspect_ratio = :equal,
    size   = (800, 700),
    legend = false,
)
savefig(plt_mask, "examples/tritronquee_mask.png")
println("  Wrote examples/tritronquee_mask.png")

println("\nDone. View the PNGs from your host filesystem.")
