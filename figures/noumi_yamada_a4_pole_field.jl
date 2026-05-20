# figures/noumi_yamada_a4_pole_field.jl
#
# The A_4^(1) Noumi–Yamada pole-field figure — the first *visual*
# demonstration that the v0.2 vector stack maps the pole field of a
# higher-order Painlevé system.  Satisfies the v0.2 PRD acceptance
# item (`docs/v0p2_plan.md` row V8a):
#
#   > An A_n^(1) pole-field figure, n ≥ 4, reproducible from
#   > `figures/`.
#
# Here `n` is the *subscript* of `A_n^(1)`: `A_4^(1)` is the `n ≥ 4`
# case, built in `NoumiYamada` as the even-parity system with index
# `n = 2` (`l = 2n = 4`, five components `f_0,…,f_4`).
#
# ## What A_4^(1) is
#
# The Noumi–Yamada `A_4^(1)` system (Noumi & Yamada 1998;
# `docs/v0p2_pillarB_noumi_yamada_findings.md` §1) is a first-order
# system for five unknowns `f_0,…,f_4` with five parameters
# `α_0,…,α_4`:
#
#     f_j' = f_j ( f_{j+1} − f_{j+2} + f_{j+3} − f_{j+4} ) + α_j ,
#
# indices mod 5.  It is one rung of the higher-order Painlevé
# hierarchy: `A_2^(1)` reduces to scalar PIV, and `A_4^(1)` is the
# first genuinely *higher*-order member — the v0.2 north-star target.
# For the `k = 1` normalisation (`Σα_j = 1`) the components satisfy
# the flow invariant `Σf_j = t` (pillar B §1.2).
#
# ## Why the shared-Q roots ARE the pole field
#
# The v0.2 vector integrator fits all five components with a *single*
# shared denominator `Q` (ADR-0019; type-II Hermite–Padé).  A vector
# meromorphic solution's components all blow up *together* at the
# same movable pole, so the roots of the one shared `Q` are *every*
# component's poles simultaneously — the pole field of the whole
# system, one consistent estimate (`src/VectorPoleField.jl`).
# `extract_poles_shared_q` roots the per-node `Q`'s, maps the roots
# to the `t`-plane, and cross-node-clusters; the result is the pole
# cloud this figure scatters.
#
# ## The IC / α / region chosen
#
# The *trivial* rational solutions of A_4^(1) (pillar B §5.2 — Types
# A/B/C: `f_j = t/5` for Type C, etc.) have **no finite poles** — a
# pole-free degeneracy that would render an empty figure.  A
# *generic* (non-rational) solution is what carries an interesting
# pole field.  We therefore use:
#
#   * α = (0.30, 0.10, 0.25, 0.15, 0.20) — five *distinct* generic
#     parameters summing to exactly 1.  Real ⇒ the A_4^(1) RHS is
#     conjugate-symmetric.
#   * f0 = (0.7, -0.3, 0.5, -0.55, -0.35) — an **all-nonzero**
#     initial condition.  The bead's V5b warning: a component
#     identically zero triggers the shared-Q degeneracy, so every
#     component is non-zero; the fifth is pinned so Σf0 = 0 exactly
#     (the constraint Σf_j = t at t0 = 0).
#   * region: the square [-2.5, 2.5]² of the complex t-plane, a
#     0.5-spaced raster — the largest region the wedge walk
#     navigates cleanly while still showing a non-trivial pole cloud.
#
# These were arrived at empirically; the all-nonzero IC is the
# load-bearing choice (it dodges the V5b degeneracy).  Every
# constant lives in `_noumi_yamada_a4_helpers.jl`, shared verbatim
# with `test/noumi_yamada_a4_figure_test.jl`.
#
# ## The figure — three panels
#
#   * Panel A — the pole field: the extracted poles as black dots in
#     the complex t-plane, over the [-3.5, 3.5]² window.  This is the
#     PRD's "A_n^(1) pole-field figure".  A faint box marks the
#     [-2.5, 2.5]² target region the walk sampled.
#
#   * Panel B — the Stage-1 visited tree: the `visited_parent` edges
#     of the wedge walk (the FW 2011 Fig 3.2 path-tree, lifted to the
#     vector setting), rooted at the IC `t0 = 0`.  Shows how the walk
#     threaded the region; step density rises near pole-dense zones.
#
#   * Panel C — ‖f(t)‖ at the visited nodes: a scatter coloured by
#     log10‖f‖, the vector-norm proxy for "distance to the next
#     pole" (the min-‖y‖ wedge criterion, `src/VectorPathNetwork.jl`).
#     Bright nodes sit near poles — visually corroborating Panel A.
#
# ## Acceptance
#
# This is a *visual* pole-field figure; there is no published
# A_4^(1) pole-field reference to pin pixel-for-pixel.  The
# quantitative obligation lives in `test/noumi_yamada_a4_figure_
# test.jl`, which reproduces this script's exact computation (via the
# shared helper) and asserts genuine invariants — the `Σf_j = t`
# constraint along the walk, a finite non-empty pole field, pole
# genuineness (`Q` vanishes there), determinism, and the equation's
# exact conjugate symmetry.
#
# References:
#   * `docs/v0p2_pillarB_noumi_yamada_findings.md` §1, §5.2.
#   * `docs/v0p2_plan.md` row V8a — the PRD `n ≥ 4` acceptance item.
#   * `src/NoumiYamada.jl`, `src/VectorPathNetwork.jl`,
#     `src/VectorPoleField.jl`.

using CairoMakie
using Printf

include(joinpath(@__DIR__, "_noumi_yamada_a4_helpers.jl"))

const OUTPNG = joinpath(@__DIR__, "output", "noumi_yamada_a4_pole_field.png")

# ----------------------------------------------------------------------
# Compute — the exact pipeline the acceptance test reproduces.
# ----------------------------------------------------------------------
@printf("A_4^(1) Noumi–Yamada pole field (v0.2 plan row V8a)\n")
@printf("  α  = %s\n", string(round.(ComplexF64.(NY_ALPHA); digits = 3)))
@printf("  f0 = %s   (Σf0 = %s)\n",
        string(round.(NY_F0; digits = 3)), string(round(sum(NY_F0); digits = 12)))
@printf("  region: [-%.1f, %.1f]^2 of the complex t-plane, step %.2f\n",
        NY_REGION_R, NY_REGION_R, NY_REGION_STEP)

t_walk = time()
sol = a4_walk()
@printf("  Stage-1 walk: %d visited nodes in %.2f s\n",
        length(sol.visited_z), time() - t_walk)

poles = a4_pole_field(sol)
@printf("  pole field: %d poles extracted (shared-Q roots)\n", length(poles))

resid = a4_constraint_residual(sol)
@printf("  max |Σf_j − t| along the walk: %.3e\n", resid)

# ‖f‖ at every visited node — the min-‖y‖ wedge proxy.
node_norm = [sqrt(sum(abs2, y)) for y in sol.visited_y]

# ----------------------------------------------------------------------
# Render — three panels.
# ----------------------------------------------------------------------
const WIN = 3.5    # plotting half-window in the t-plane

fig = Figure(size = (1500, 540))
Label(fig[0, 1:3],
      "A_4⁽¹⁾ Noumi–Yamada pole field — generic α = (0.30, 0.10, 0.25, 0.15, 0.20), " *
      "all-nonzero IC, complex-t region [−2.5, 2.5]²";
      fontsize = 15, padding = (0, 0, 8, 8))

# ---- Panel A: the pole field ----------------------------------------
axA = Axis(fig[1, 1];
           xlabel = "Re t", ylabel = "Im t",
           title  = "A. Pole field — roots of the shared denominator Q",
           subtitle = "$(length(poles)) poles; PRD n ≥ 4 acceptance",
           titlesize = 13, subtitlesize = 10,
           aspect = DataAspect(),
           limits = (-WIN, WIN, -WIN, WIN))
# faint axes + the [-R, R]^2 target box
lines!(axA, [-WIN, WIN], [0.0, 0.0]; color = (:gray, 0.35), linewidth = 0.5)
lines!(axA, [0.0, 0.0], [-WIN, WIN]; color = (:gray, 0.35), linewidth = 0.5)
let R = NY_REGION_R
    lines!(axA, [-R, R, R, -R, -R], [-R, -R, R, R, -R];
           color = (:steelblue, 0.6), linewidth = 0.8, linestyle = :dash)
end
scatter!(axA, Float64.(real.(poles)), Float64.(imag.(poles));
         color = :black, markersize = 7)
# IC marker at t0 = 0
scatter!(axA, [Float64(real(NY_T0))], [Float64(imag(NY_T0))];
         color = :red, markersize = 9, marker = :star5)

# ---- Panel B: the Stage-1 visited tree ------------------------------
axB = Axis(fig[1, 2];
           xlabel = "Re t", ylabel = "Im t",
           title  = "B. Stage-1 visited tree (visited_parent edges)",
           subtitle = "$(length(sol.visited_z)) nodes; wedge walk from t₀ = 0",
           titlesize = 13, subtitlesize = 10,
           aspect = DataAspect(),
           limits = (-WIN, WIN, -WIN, WIN))
seg = Point2f[]
for k in 2:length(sol.visited_z)
    p = sol.visited_parent[k]
    p == 0 && continue
    push!(seg, Point2f(real(sol.visited_z[p]), imag(sol.visited_z[p])))
    push!(seg, Point2f(real(sol.visited_z[k]), imag(sol.visited_z[k])))
end
linesegments!(axB, seg; color = (:black, 0.6), linewidth = 0.5)
scatter!(axB, Float64.(real.(sol.visited_z)), Float64.(imag.(sol.visited_z));
         color = :black, markersize = 2.0)
scatter!(axB, [Float64(real(NY_T0))], [Float64(imag(NY_T0))];
         color = :red, markersize = 9, marker = :star5)

# ---- Panel C: ‖f(t)‖ at visited nodes -------------------------------
axC = Axis(fig[1, 3];
           xlabel = "Re t", ylabel = "Im t",
           title  = "C. log₁₀‖f(t)‖ at visited nodes",
           subtitle = "min-‖y‖ wedge proxy: bright ⇒ near a pole",
           titlesize = 13, subtitlesize = 10,
           aspect = DataAspect(),
           limits = (-WIN, WIN, -WIN, WIN))
lines!(axC, [-WIN, WIN], [0.0, 0.0]; color = (:gray, 0.35), linewidth = 0.5)
lines!(axC, [0.0, 0.0], [-WIN, WIN]; color = (:gray, 0.35), linewidth = 0.5)
scC = scatter!(axC,
               Float64.(real.(sol.visited_z)), Float64.(imag.(sol.visited_z));
               color = log10.(node_norm), colormap = :viridis, markersize = 6)
Colorbar(fig[1, 4], scC; label = "log₁₀‖f‖")
# overlay the extracted poles for visual cross-check
scatter!(axC, Float64.(real.(poles)), Float64.(imag.(poles));
         color = :red, markersize = 4, marker = :xcross)

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
