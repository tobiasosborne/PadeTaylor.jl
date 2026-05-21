# figures/kkg_pi2_tritronquee_pole_field.jl
#
# The P_I^(2) tritronquée + pole-field figure — the final *visual*
# deliverable of the v0.2 vector epic (bead `padetaylor-0ln.18`,
# v0.2 plan row V8b).  It reproduces the structure of the
# Kapaev–Klein–Grava (KKG) 2015 tritronquée study, Figs. 7.4/7.5
# (`references/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41.pdf`):
# the distinguished tritronquée solution `V_0` of the fourth-order
# Painlevé-I-hierarchy equation P_I^(2), and the pole field it carries
# in the wedge straddling the positive real axis.
#
# ## The equation, the solution, the (sign-corrected) IC
#
# P_I^(2) is the fourth-order member of the Painlevé-I hierarchy
# (KKG 2015 eq. (1.1); pillar C §1,
# `docs/v0p2_pillarC_painleve_hierarchy_findings.md`), here at `t = 0`:
#
#     u'''' + 10 u'^2 + 20 u u'' + 40 (u^3 - 6 t u + 6 x) = 0.
#
# Its **tritronquée** `V_0` (KKG 2015 §7) is the solution pole-free
# over a sector of angular width ≈ 270° of the complex `x`-plane
# (Claeys 2011), with a pole field only in the remaining wedge
# straddling `arg x ≈ 0` (the positive real axis; pillar C items 1–2).
# On the negative real axis it is real, positive and monotone, with
# the algebraic asymptotics `u → +(6|x|)^{1/3}` as `x → -∞` (KKG 2015
# eq. (1.3); pillar C item 4).
#
# The initial condition is the KKG asymptotic seed
# `pI2_tritronquee_ic(x0; n_terms = 2)` — the sign-corrected library
# seed.  The sign is decisive: the dominant balance `40(u^3+6x)=0`
# gives `u^3 = 6|x| > 0` on `x < 0`, so the real branch is POSITIVE
# (`u(-20) ≈ +4.93`).  Seeding the negative branch is not a P_I^(2)
# solution at all (ODE residual `≈ -480|x|`).
#
# ## How the figure is computed — two stages (shared kernel)
#
# Every load-bearing constant and computation lives in the Makie-free
# kernel `figures/_kkg_pi2_helpers.jl`, shared verbatim with the
# acceptance test `test/kkg_pi2_figure_test.jl`.  This script only
# renders.
#
# **Stage A (Panel A) — the tritronquée via the vector BVP.**  Forward
# IVP integration of a tritronquée diverges (it is a separatrix), so a
# global Chebyshev-collocation BVP is mandatory (`src/VectorBVP.jl`).
# The companion 4-vector P_I^(2) system is collocated on a
# negative-real-axis segment `[-20, -2]` with a 2+2 split boundary
# condition pinning `(u, u')` at both endpoints from the asymptotic
# seed.  This is the recipe verified by probe `padetaylor-0ln.29`
# (`external/probes/pi2-tritronquee-bvp/RECIPE.md`): `vector_bvp_solve`
# converges in 2–3 Newton iterations, ODE residual `~1.5e-11`, matching
# the KKG `n_terms = 2` truncated series to `~2.1e-4` on `[-20,-2]`.
#
# **Stage B (Panel B) — the pole-field march.**  A vector path-network
# wedge walk is seeded from a BVP-derived state near the right end of
# the negative-real segment and marched toward a fan of targets
# covering the pole-rich wedge near `arg x ≈ 0`;
# `extract_poles_shared_q` then roots the per-node shared denominators
# `Q` (one shared `Q` for all four companion components — ADR-0019) and
# cross-node-clusters.  The result is the pole scatter of Panel B.
#
# ## V1 corners (CLAUDE.md Rule 9 — documented, with the forcing
# condition)
#
#   * **The wedge step `h`.**  The tritronquée's pole field is dense in
#     the wedge; a coarse step overshoots a pole onto a blown-up state
#     where the per-node canonical-`Q` store cannot be built.  `h` is
#     therefore hand-tuned (`h = 0.1` — see `_kkg_pi2_helpers.jl`).  The
#     principled fix — a shared-`Q`-root-distance wedge criterion — is
#     deferred bead `padetaylor-0ln.23`.
#
#   * **The full KKG Fig 7.4/7.5 Re/Im surface.**  KKG plot `Re V_0`,
#     `Im V_0` as filled surfaces over `[-20,20]^2`.  This figure
#     scatters the pole *locations* only — a dense filled heatmap needs
#     the Stage-2 fine-grid fill, the deferred bead `padetaylor-0ln.28`.
#
# ## Acceptance
#
# Panel A is the quantitative deliverable — the BVP recipe is proven
# (the probe RECIPE.md numbers).  Panel B's pole-field march is the
# V8b open risk; the kernel `kkg_pole_field` is written to report a
# march failure honestly (`march_ok = false`) rather than fabricate
# poles.  The quantitative obligations are pinned in
# `test/kkg_pi2_figure_test.jl` (`PI2F.*`): BVP convergence, the
# tritronquée sign+magnitude (the regression guard for the sign bug),
# the P_I^(2) ODE residual at interior nodes, a non-empty pole field
# in the wedge, the pole-free sector, and determinism.
#
# References:
#   * `references/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41.pdf`
#     — KKG 2015: eq. (1.1), §7, Figs. 7.4–7.5.
#   * `external/probes/pi2-tritronquee-bvp/RECIPE.md` — the verified
#     2+2 (u,u') BVP recipe + numbers.
#   * `docs/v0p2_pillarC_painleve_hierarchy_findings.md` — the equation,
#     the sign-corrected IC, the pole-free sector.
#   * `docs/v0p2_plan.md` row V8b — the PRD acceptance item.
#   * `src/VectorBVP.jl`, `src/VectorPathNetwork.jl`,
#     `src/VectorPoleField.jl`, `src/PainleveHierarchy.jl`.

using CairoMakie
using Printf

include(joinpath(@__DIR__, "_kkg_pi2_helpers.jl"))

const OUTPNG = joinpath(@__DIR__, "output",
                        "kkg_pi2_tritronquee_pole_field.png")

# ----------------------------------------------------------------------
# Compute — the exact pipeline the acceptance test reproduces.
# ----------------------------------------------------------------------
@printf("P_I^(2) tritronquée + pole field (v0.2 plan row V8b; KKG 2015 §7)\n")
@printf("  BVP segment: [%.1f, %.1f] on the negative real axis, N = %d\n",
        real(KKG_X_L), real(KKG_X_R), KKG_BVP_N)

t_bvp = time()
res = kkg_pi2_compute()
@printf("  BVP solve: %d Newton iterations, residual_inf = %.3e, %.2f s\n",
        res.bvp.iterations, res.bvp.residual_inf, time() - t_bvp)

# Quantitative check against the KKG n_terms=2 asymptotic series.
kkg_err = maximum(abs(u - kkg_asymptotic(x))
                  for (x, u) in zip(res.xs, res.us))
@printf("  tritronquée: u(-15) = %.6f  (KKG (6·15)^(1/3) = %.6f)\n",
        res.us[findmin(abs.(res.xs .- (-15.0)))[2]], cbrt(90.0))
@printf("  max |u_BVP − u_KKG(n=2)| over [%.0f,%.0f] = %.3e\n",
        real(KKG_X_L), real(KKG_X_R), kkg_err)

@printf("  Stage-B pole-field march: %s\n", res.message)

# ----------------------------------------------------------------------
# Render — two panels.
# ----------------------------------------------------------------------
fig = Figure(size = (1180, 520))
Label(fig[0, 1:2],
      "P_I⁽²⁾ tritronquée V₀ (t = 0) — BVP solution on the negative real " *
      "axis + pole-field march into the arg x ≈ 0 wedge  [KKG 2015 Figs 7.4/7.5]";
      fontsize = 14, padding = (0, 0, 8, 8))

# ---- Panel A: the BVP-computed tritronquée on the negative real axis -
axA = Axis(fig[1, 1];
           xlabel = "x", ylabel = "Re u(x)",
           title  = "A. Tritronquée on the negative real axis (vector BVP)",
           subtitle = @sprintf("2+2 (u,u') BVP, N = %d, residual %.1e — u > 0, monotone",
                                KKG_BVP_N, res.bvp.residual_inf),
           titlesize = 13, subtitlesize = 10)
# The BVP solution.
lines!(axA, res.xs, res.us; color = :navy, linewidth = 2.2,
       label = "u_BVP(x)")
# The KKG n_terms=2 asymptotic series, for visual cross-check.
lines!(axA, res.xs, [kkg_asymptotic(x) for x in res.xs];
       color = :darkorange, linewidth = 1.2, linestyle = :dash,
       label = "KKG asymptotic  Y + Y⁻⁶")
axislegend(axA; position = :lt, framevisible = false, labelsize = 10)

# ---- Panel B: the pole-field scatter --------------------------------
axB = Axis(fig[1, 2];
           xlabel = "Re x", ylabel = "Im x",
           title  = "B. Pole field — roots of the shared denominator Q",
           subtitle = res.march_ok ?
               @sprintf("%d poles in the arg x ≈ 0 wedge; %d path-network nodes",
                        length(res.poles), length(res.walk.visited_z)) :
               "Stage-B march BLOCKED — see console / figure note",
           titlesize = 13, subtitlesize = 10,
           aspect = DataAspect())
# faint axes
let W = 10.0
    lines!(axB, [-W, W], [0.0, 0.0]; color = (:gray, 0.35), linewidth = 0.5)
    lines!(axB, [0.0, 0.0], [-W, W]; color = (:gray, 0.35), linewidth = 0.5)
end
if res.march_ok
    # path-network visited nodes, faint behind
    scatter!(axB,
             Float64.(real.(res.walk.visited_z)),
             Float64.(imag.(res.walk.visited_z));
             color = (:steelblue, 0.35), markersize = 3.0,
             label = "path-network nodes")
    # the extracted pole field
    scatter!(axB,
             Float64.(real.(res.poles)), Float64.(imag.(res.poles));
             color = :crimson, markersize = 9, marker = :circle,
             label = "extracted poles")
    axislegend(axB; position = :lt, framevisible = false, labelsize = 9)
else
    text!(axB, 0.0, 0.0;
          text = "pole-field march blocked\n(" * res.message * ")",
          align = (:center, :center), fontsize = 11, color = :crimson)
end
# the BVP seed point of the march
scatter!(axB, [Float64(real(KKG_Z_SEED))], [Float64(imag(KKG_Z_SEED))];
         color = :black, markersize = 10, marker = :star5)

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
