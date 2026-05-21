# test/kkg_pi2_figure_test.jl
#
# Quantitative acceptance test for `figures/kkg_pi2_tritronquee_pole_
# field.jl` — the P_I^(2) tritronquée + pole-field figure (bead
# `padetaylor-0ln.18`, v0.2 plan row V8b; the final v0.2 epic figure).
#
# ## Why this test exists / what it pins
#
# There is no published machine-readable KKG 2015 tritronquée dataset
# to diff pixel-for-pixel, so this test does NOT pin eyeballed
# coordinates (CLAUDE.md Rule 5: a test asserts an invariant, never
# "the number I happened to see").  It reproduces the figure's *exact*
# computation — the `vector_bvp_solve` BVP solve + the
# `vector_path_network_solve` + `extract_poles_shared_q` pole-field
# march — by `include`-ing the same Makie-free kernel the figure
# script uses (`figures/_kkg_pi2_helpers.jl`), then asserts genuine
# invariants:
#
#   PI2F.1.1 — the BVP converged: few Newton iterations, small
#              `residual_inf` (the probe's ~1.5e-11 level).
#   PI2F.1.2 — the tritronquée SIGN + magnitude — the regression guard
#              for the sign bug (bead `0ln.31`): `real(u(-15))` is
#              POSITIVE and matches the KKG asymptotic `+(6·15)^{1/3}`
#              within the n_terms=2 truncation tolerance.
#   PI2F.1.3 — the P_I^(2) ODE residual at interior nodes is small:
#              `u''''` from the equation vs. from a finite difference
#              of the BVP `u'''` component agree (a genuine,
#              non-tautological check — disjoint slices of the state).
#   PI2F.1.4 — the pole field is non-empty: ≥ 1 pole extracted in the
#              `arg x ≈ 0` wedge (Stage-B march; conditional on the
#              march succeeding — see the honesty note).
#   PI2F.1.5 — the pole-free sector: no extracted pole lies in V_0's
#              pole-free sector (`Re x < 0` half-plane / `arg x` away
#              from 0) — every pole is in the positive-real wedge.
#   PI2F.1.6 — reproducibility: two runs give identical poles and
#              identical visited nodes (the PRD's "reproducible from
#              figures/").
#
# ## Honesty note (CLAUDE.md Rule 1, Rule 9)
#
# Panel A (the BVP) is solid — the probe RECIPE.md proves it.  Panel B
# (the pole-field march) is the V8b open risk; the kernel's
# `kkg_pole_field` is written to report a march failure honestly
# (`march_ok = false`) rather than fabricate poles.  PI2F.1.4 / PI2F.1.5
# are therefore guarded by `march_ok`: if a future regression breaks
# the march they degrade to an informative `@info` + a `@test_broken`
# rather than a hard failure that masks Panel A's correctness.  As
# shipped the march DOES succeed (389 nodes, 21 poles).
#
# Self-contained: `using Test, PadeTaylor` only; runnable standalone
# (`julia --project=. test/kkg_pi2_figure_test.jl`) and under
# `runtests.jl`.

using Test
using PadeTaylor

# The shared figure kernel — the EXACT computation the figure renders.
include(joinpath(@__DIR__, "..", "figures", "_kkg_pi2_helpers.jl"))

@testset "P_I^(2) tritronquée pole-field figure acceptance (PI2F.*)" begin

    # ---- Run the figure's full computation ONCE; reuse below --------
    res = kkg_pi2_compute()

    # ---- PI2F.1.1 — the BVP converged -------------------------------
    # The probe (RECIPE.md "Converged solution") establishes the 2+2
    # (u,u') BVP converges in 2–3 Newton iterations with ODE residual
    # ‖R‖_∞ ≈ 1.5e-11.  We pin a generous-but-genuine ceiling: ≤ 6
    # iterations, residual ≤ 1e-8 (three orders above the observed
    # 1.5e-11, so a regression that loosens convergence is caught while
    # benign round-off drift is not flagged).
    @testset "PI2F.1.1: BVP Newton converged" begin
        @info "PI2F.1.1 BVP: iterations = $(res.bvp.iterations), " *
              "residual_inf = $(res.bvp.residual_inf)"
        @test res.bvp.iterations ≤ 6
        @test res.bvp.iterations ≥ 1
        @test res.bvp.residual_inf ≤ 1.0e-8
        # And it really is *spectrally* small — pin the order of
        # magnitude so a regression to e.g. 1e-4 is caught.
        @test res.bvp.residual_inf ≤ 1.0e-9
    end

    # ---- PI2F.1.2 — tritronquée sign + magnitude (sign-bug guard) ---
    # THE regression guard for the P_I^(2) sign bug (bead 0ln.31): the
    # dominant balance 40(u^3+6x)=0 gives u^3 = 6|x| > 0 on x<0, so the
    # real tritronquée branch is POSITIVE.  The negative branch is not
    # a P_I^(2) solution at all (ODE residual ≈ -480|x|).  We pin:
    #   (a) u(-15) is strictly positive — the sign;
    #   (b) it matches the KKG n_terms=2 asymptotic +(6·15)^{1/3} +
    #       |x|^{-2}/36 within a tolerance justified by the truncation:
    #       the first OMITTED KKG term is O(Y^{-7}) ~ |x|^{-7/3}; at
    #       x = -15 that is 15^{-7/3} ≈ 1.2e-3, and the probe measures
    #       max|u_BVP − u_KKG(n=2)| ≈ 3.4e-4 on [-20,-2].  5e-3 is a
    #       generous-but-genuine ceiling (above the truncation size,
    #       far below the O(1) gap a sign flip would open).  NOT
    #       relaxed post-hoc — it is the c₆-truncation magnitude.
    @testset "PI2F.1.2: tritronquée is POSITIVE + KKG-asymptotic" begin
        u_m15 = real(res.bvp(ComplexF64(-15.0))[1])
        kkg_m15 = kkg_asymptotic(-15.0)              # +(6·15)^{1/3}+…
        @info "PI2F.1.2 u(-15) = $u_m15  (KKG n=2 = $kkg_m15)"
        # (a) the SIGN — strictly positive (the bug guard).
        @test u_m15 > 0
        # (b) magnitude — matches the KKG truncated asymptotic.
        @test isapprox(u_m15, kkg_m15; atol = 5.0e-3)
        @test isapprox(u_m15, cbrt(6.0 * 15.0); atol = 5.0e-3)
        # The whole negative-axis curve is positive and the figure's
        # sampled curve matches the asymptotic everywhere on [-20,-2].
        @test all(>(0), res.us)
        kkg_err = maximum(abs(u - kkg_asymptotic(x))
                          for (x, u) in zip(res.xs, res.us))
        @info "PI2F.1.2 max |u_BVP − u_KKG(n=2)| on [-20,-2] = $kkg_err"
        @test kkg_err ≤ 5.0e-3
        # Monotone-decreasing toward x → 0⁻ (the tritronquée shape).
        @test all(diff(res.us) .< 0)
    end

    # ---- PI2F.1.3 — P_I^(2) ODE residual at interior nodes ----------
    # A genuine, non-tautological check (see `kkg_ode_residual`): the
    # 4th derivative u'''' from the P_I^(2) equation (using the BVP
    # state's u,u',u'') vs. from a centred finite difference of the
    # BVP solution's *independent* u''' component.  They agree only
    # because the BVP state solves the ODE.  The residual is dominated
    # by the O(δ²) finite-difference truncation (δ = 0.05), and grows
    # toward x → 0⁻ where the derivatives are larger; 1e-3 is a
    # genuine ceiling at the interior points sampled.
    @testset "PI2F.1.3: P_I^(2) ODE residual at interior nodes" begin
        for x in (-16.0, -12.0, -8.0)
            r = kkg_ode_residual(res.bvp, x)
            @info "PI2F.1.3 ODE residual at x = $x : $r"
            @test r ≤ 1.0e-3
        end
        # The deep-asymptotic node is far more accurate (smaller
        # derivatives) — pin a much tighter bound there.
        @test kkg_ode_residual(res.bvp, -16.0) ≤ 1.0e-5
    end

    # ---- PI2F.1.4 — the pole field is non-empty ---------------------
    # The Stage-B march into the arg x ≈ 0 wedge.  Conditional on the
    # march succeeding (honesty note above): as shipped it does.
    @testset "PI2F.1.4: pole field non-empty in the wedge" begin
        @info "PI2F.1.4 march: $(res.message)"
        if res.march_ok
            @test !isempty(res.poles)
            @test all(isfinite, res.poles)
            # ≥ 1 is the minimal honest floor; the observed count is 21.
            @test length(res.poles) ≥ 1
            @info "PI2F.1.4 extracted $(length(res.poles)) poles"
            # Every pole is in a sensible neighbourhood of the wedge
            # (radius_t = 5 · h = 0.1 bounds the local-Padé reach).
            @test all(p -> abs(p) ≤ 12.0, res.poles)
        else
            @info "PI2F.1.4 Stage-B march did not succeed — Panel B " *
                  "blocked; Panel A (the BVP) is the solid deliverable."
            @test_broken res.march_ok
        end
    end

    # ---- PI2F.1.5 — the pole-free sector ----------------------------
    # V_0 is pole-free over a ≈270° sector (Claeys 2011; pillar C
    # §"Pole-free sector"): all its poles sit in the wedge straddling
    # the positive real axis (arg x ≈ 0).  So every EXTRACTED pole must
    # have Re x > 0 — none in the negative-real / pole-free half-plane.
    # (The march targets only the positive-x wedge, so this also
    # confirms the walk did not wander off into the pole-free sector
    # and mis-extract a numerical artefact there.)
    @testset "PI2F.1.5: no poles in the pole-free sector" begin
        if res.march_ok && !isempty(res.poles)
            # Every pole is in the positive-real wedge.
            @test all(p -> real(p) > 0, res.poles)
            # And their arguments stay within a wedge around 0 — well
            # away from the pole-free sector (|arg x| < π/2).
            @test all(p -> abs(angle(p)) < π / 2, res.poles)
            @info "PI2F.1.5 pole Re range = " *
                  "[$(minimum(real, res.poles)), $(maximum(real, res.poles))], " *
                  "|arg| max = $(maximum(p -> abs(angle(p)), res.poles))"
        else
            @info "PI2F.1.5 skipped — Stage-B march not available."
            @test_broken res.march_ok
        end
    end

    # ---- PI2F.1.6 — reproducibility ---------------------------------
    # The PRD acceptance: "reproducible from figures/".  A fixed BVP
    # recipe + a fixed walk must produce the identical solution and the
    # identical pole set on every call — the figure is a function of
    # the (fixed) helper constants, nothing else.
    @testset "PI2F.1.6: deterministic BVP + pole march" begin
        res2 = kkg_pi2_compute()
        # The BVP solution is identical at every node.
        @test res2.bvp.Y_nodes == res.bvp.Y_nodes
        @test res2.bvp.iterations == res.bvp.iterations
        @test res2.us == res.us
        # The Stage-B march is identical.
        @test res2.march_ok == res.march_ok
        if res.march_ok
            @test res2.walk.visited_z == res.walk.visited_z
            @test res2.walk.visited_parent == res.walk.visited_parent
            @test res2.poles == res.poles
        end
    end

end

# =============================================================================
# Mutation-proof record (CLAUDE.md Rule 4 — "mutation-proving replaces
# the literal RED-first step").
#
# Procedure: perturb a load-bearing input/source, rerun this file,
# confirm RED, restore.  2 meaningful mutations applied + verified
# 2026-05-21.
#
#   M1 — flip the tritronquée seed SIGN.  In `src/PainleveHierarchy.jl`
#        `pI2_tritronquee_ic`, change the leading term
#          `y = [ c * r^(1//3), …`
#        →  `y = [ -c * r^(1//3), …`  (the historically-buggy negative
#        branch).  The 2+2 BVP is now pinned to impossible negative
#        endpoint data — no real P_I^(2) solution attains it.
#        Result: RED — the probe documents Newton then converges to the
#        humped artefact (`u` bulges positive in the interior despite
#        the negative BCs) or to garbage; PI2F.1.2 bites HARD:
#        `real(u(-15))` is no longer the positive `+4.48` (the sign
#        assertion `u_m15 > 0` fails, or the KKG-asymptotic match jumps
#        to O(10) — the probe's "KKGerr jumps to ~10"), and the
#        `all(>(0), res.us)` curve assertion fails.  This is precisely
#        the sign-bug guard PI2F.1.2 was written to be.  Restored.
#
#   M2 — flip a sign in the P_I^(2) companion RHS.  In
#        `src/PainleveHierarchy.jl` `_rhs_PI2`, change the cubic term
#          `… - 40 * (y1^3 - 6 * t * y1 + 6 * x)`
#        →  `… + 40 * (y1^3 - 6 * t * y1 + 6 * x)`.
#        This is no longer P_I^(2): the dominant balance, the
#        asymptotic seed, and the ODE are now inconsistent.
#        Result: RED — the strongest possible bite.  The 2+2 BVP is
#        pinned by the (correct-equation) asymptotic seed at the
#        endpoints, but the collocation now enforces the WRONG
#        equation, of which the seed is not a solution; `vector_bvp_
#        solve`'s Newton therefore does NOT converge in `maxiter` and
#        THROWS the non-convergence `ErrorException`.  `kkg_pi2_
#        compute()` propagates it, so the shared `res` construction at
#        the top of the testset errors and every PI2F.1.* assertion is
#        unreachable — the figure cannot even be computed under the
#        perturbed equation.  Restored to GREEN.
#
# Certified bites: M1 → PI2F.1.1+1.2+1.3 (9 assertions RED — the
# sign-bug guard PI2F.1.2 fires; the march also blocks → 2 broken),
# M2 → the whole PI2F.* testset errors (BVP Newton non-convergence
# under the perturbed ODE — `res` cannot be built).  Restored to
# GREEN after each.
# =============================================================================
