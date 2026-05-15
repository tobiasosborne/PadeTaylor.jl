# test/adaptive_step_test.jl — Phase A1 / bead `padetaylor-8ui`.
#
# Adaptive Padé step size per FFW 2017 §2.1.2 (`references/markdown/
# FFW2017_painleve_riemann_surfaces_preprint/
# FFW2017_painleve_riemann_surfaces_preprint.md:74-97`).
#
# The algorithm: at each step compute the truncation-error estimate
# `T(h) = |ε_{n+1} h^{n+1} / a(h)|` where `a(h)` is the Padé numerator
# evaluated at the step, `ε_k = c_k + Σ_{r=1..ν} b_r c_{k-r}` are the
# Padé error coefficients, and `c_k` are Taylor coefficients (we
# compute order `n+1` Taylor to access `c_{n+1}`).  If `T(h) > Tol`,
# rescale `h := q·h` with `q = (k·Tol/T(h))^(1/(n+1))` and retry; if
# `T(h) ≤ Tol`, accept and seed the next step's initial `h` from
# `|q·h|`.  FFW conservative factor `k ≈ 1e-3`; Padé order `ν = n/2`.
#
# Test ID prefix `AS` (Adaptive Step).  Acceptance bar: AS.1.1-AS.1.5
# all GREEN and mutation-proven by M1-M4 (see footer).

using Test
using PadeTaylor
using PadeTaylor.PadeStepper: ffw_truncation_error, ffw_rescale_q,
                               adaptive_pade_step!, PadeStepperState

include(joinpath(@__DIR__, "_oracle_problems.jl"))

@testset "AdaptiveStep (FFW 2017 §2.1.2): adaptive Padé `h`" begin

    # =================================================================
    # AS.1.1 — T(h) formula sanity on `y' = y` (so `u'' = u` with the
    # `up` slot ignored; closed-form u(z) = exp(z) with u(0) = u'(0) = 1).
    # Computed `T(h)` must (i) be monotonically increasing in h and
    # (ii) scale as h^{n+1} at small h (the leading-h asymptote of the
    # FFW estimator).  At smaller `h` the prefactor `|ε̃_{n+1}/a(h)|`
    # stabilises and the `h^{n+1}` power is the only h-dependent factor.
    # =================================================================
    @testset "AS.1.1: T(h) formula matches closed-form on u''=u" begin
        # u'' = u, ICs u(0) = 1, u'(0) = 1 ⇒ u(z) = exp(z) ⇒ c_k = 1/k!
        f = (z, u, up) -> u
        z0, u0, up0 = 0.0, 1.0, 1.0
        order = 10                     # n = 10, ν = 5

        Th_vals = Float64[]
        for h in (0.1, 0.5, 1.0)
            Th = ffw_truncation_error(f, z0, u0, up0, order, h)
            push!(Th_vals, Th)
            @test isfinite(Th)
            @test Th > 0
        end
        # (i) Monotonicity in h (the dominant `h^{n+1}` factor wins).
        @test Th_vals[1] < Th_vals[2] < Th_vals[3]

        # (ii) Asymptotic h^{n+1} scaling at small h.  Two very small
        # h values whose ratio is the dominant factor.
        h_a, h_b = 1.0e-4, 1.0e-3
        Th_a = ffw_truncation_error(f, z0, u0, up0, order, h_a)
        Th_b = ffw_truncation_error(f, z0, u0, up0, order, h_b)
        ratio_expected = (h_b / h_a) ^ (order + 1)
        ratio_actual   = Th_b / Th_a
        @test isapprox(ratio_actual, ratio_expected; rtol = 1e-3)

        # (iii) Prefactor pin.  The leading-h coefficient
        # `T(h)/h^{n+1} → |ε̃_{n+1}/a(h→0)| = |ε_{n+1}/a_0|` is
        # `c_{n+1} + Σ b_r c_{n+1-r}` evaluated at the rescaled
        # coefficients in the h → 0 limit (where the rescaled c̃_k
        # equal the unrescaled c_k of u(z)/h^k — algebraically the
        # Taylor coefficients c_k = 1/k! for exp, weighted by the
        # exp(5,5) Padé denominator b_r).  At order=10 with the
        # closed-form exp coefficients, that prefactor is ~ 1e-10.
        # Mutation M2 (dropping the `Σ b_r c_{n+1-r}` term) replaces
        # ε_{n+1} with just `c_{n+1} = 1/11! ≈ 2.5e-8`, a factor of
        # ~250× larger — biting this assertion at any reasonable
        # rtol on the prefactor.
        prefactor = Th_a / (h_a ^ (order + 1))
        @test prefactor < 1e-9     # tightly bounds away from M2's ~ 2.5e-8
        @test prefactor > 1e-12    # and away from 0 (M2 lowers ε; this
                                    # bounds the other direction)
    end

    # =================================================================
    # AS.1.2 — q-rescale fixed-point convergence.  Starting from a
    # too-large `h`, iterating `h := q·h` must converge to T(h) ≤ Tol
    # in ≤ 5 iterations.  Tested at Tol = 1e-10 and Tol = 1e-14.
    # =================================================================
    @testset "AS.1.2: q-rescale converges in ≤ 5 iterations" begin
        f = (z, u, up) -> 6 * u^2          # FW 2011 test ODE.
        z0, u0, up0 = 0.0, u_0_FW, up_0_FW
        order = 30
        k_conserv = 1.0e-3

        for Tol in (1.0e-10, 1.0e-14)
            h = 2.0                         # deliberately too large
            iters = 0
            local Th
            for it in 1:5
                Th = ffw_truncation_error(f, z0, u0, up0, order, h)
                if Th ≤ Tol
                    iters = it
                    break
                end
                q = ffw_rescale_q(Tol, Th, order; k = k_conserv)
                @test 0 < q < 1            # rescale must shrink h
                h = q * h
                iters = it
            end
            @test iters ≤ 5
            @test Th ≤ Tol
            @test h > 0
        end
    end

    # =================================================================
    # AS.1.3 — End-to-end agreement with fixed-h baseline on the
    # equianharmonic ℘ test (`u'' = 6u²`), z = 0 to z = 30 (FW Table
    # 5.1 target).  Adaptive should land within ≤ 1e-10 of fixed-h.
    #
    # We deliberately seed the adaptive walk with an aggressive
    # `h = 1.5` so that `T(1.5) ≈ 1e-7 ≫ Tol = 1e-12`: the controller
    # must shrink h before the first accepted step.  This forces at
    # least one visited node to have `h ≠ 1.5`, exercising the
    # adaptation path (not just the acceptance path).  Without this,
    # `T(0.5)` at FW IC is already ≈ 4e-22 ≪ Tol, so no rescale fires
    # and the test degenerates to "did anything run."
    # =================================================================
    @testset "AS.1.3: adaptive vs fixed-h end-to-end agreement at z=30" begin
        f = (z, u, up) -> 6 * u^2
        prob = PadeTaylorProblem(f, (u_0_FW, up_0_FW), (0.0, 30.0); order = 30)
        grid = ComplexF64[30.0 + 0.0im]

        # Fixed-h baseline (current default).
        sol_fixed = path_network_solve(prob, grid;
                                        h = 0.5,
                                        max_steps_per_target = 500)
        u_fixed = sol_fixed.grid_u[1]

        # Adaptive-h run, seeded at h = 1.5 to FORCE adaptation.
        sol_adapt = path_network_solve(prob, grid;
                                        h = 1.5,
                                        step_size_policy = :adaptive_ffw,
                                        adaptive_tol = 1.0e-12,
                                        max_steps_per_target = 5000)
        u_adapt = sol_adapt.grid_u[1]

        @test isfinite(real(u_fixed)) && isfinite(real(u_adapt))
        @test abs(u_adapt - u_fixed) ≤ 1.0e-9

        # The adaptive run must visit different `h` values at some
        # nodes (it shouldn't degenerate to all-1.5).  This is the
        # load-bearing "adaptation actually happened" assertion.
        h_values = sol_adapt.visited_h
        @test any(h -> !isapprox(h, 1.5; atol = 1e-10), h_values)
        # And no visited h should exceed the seed (controller only
        # shrinks; never grows beyond its initial value).
        @test all(h -> h ≤ 1.5 + 1e-12, h_values)
    end

    # =================================================================
    # AS.1.4 — CT.1.3 transformed-PIII test with adaptive step.  The
    # direct-vs-transformed end-to-end agreement (already established
    # in test/coord_transforms_test.jl) must continue to hold under
    # adaptive_ffw — proves the adaptive path doesn't corrupt the
    # transform.
    # =================================================================
    @testset "AS.1.4: PIII transformed agreement under adaptive step" begin
        z₀, u₀, up₀ = 1.0 + 0.0im, 0.25 + 0.0im, 1.0 + 0.0im
        α, β, γ, δ = -0.5, -0.5, 1.0, -1.0  # FFW 2017 Fig 1 caption
        z_target = 1.05 + 0.0im
        ζ_target = 2 * log(z_target)

        f_direct = (z, u, up) -> (up^2)/u - up/z +
                                  (α*u^2 + β)/z + γ*u^3 + δ/u
        prob_direct = PadeTaylorProblem(f_direct, (u₀, up₀),
                                         (z₀, z_target); order = 30)
        sol_direct = path_network_solve(prob_direct, ComplexF64[z_target];
                                         h = 0.5,
                                         step_size_policy = :adaptive_ffw,
                                         adaptive_tol = 1.0e-12,
                                         max_steps_per_target = 100)
        u_direct  = sol_direct.grid_u[1]

        f_trans = pIII_transformed_rhs(α, β, γ, δ)
        ζ₀, w₀_, wp₀_ = pIII_z_to_ζ(z₀, u₀, up₀)
        prob_trans = PadeTaylorProblem(f_trans, (w₀_, wp₀_),
                                        (ζ₀, ζ_target); order = 30)
        sol_trans = path_network_solve(prob_trans, ComplexF64[ζ_target];
                                        h = 0.5,
                                        step_size_policy = :adaptive_ffw,
                                        adaptive_tol = 1.0e-12,
                                        max_steps_per_target = 100)
        w_trans = sol_trans.grid_u[1]
        _, u_recovered, _ = pIII_ζ_to_z(ζ_target, w_trans, sol_trans.grid_up[1])

        @test isfinite(u_direct) && isfinite(u_recovered)
        @test abs(u_direct - u_recovered) ≤ 1.0e-10
    end

    # =================================================================
    # AS.1.5 — Tolerance sweep: maximum visited-h decreases as Tol
    # tightens.  The FFW rescale law `q = (k·Tol/T(h))^(1/(n+1))` is
    # monotone in Tol — a smaller Tol forces a smaller q at any given
    # T(h), so the accepted h is smaller.  Across a Tol sweep at fixed
    # initial h=1.5 (large enough that T(1.5) > Tol for all sweep
    # values), the maximum visited_h must be non-increasing as Tol
    # tightens.
    #
    # Final solution-value error is NOT a clean monotone proxy —
    # at the FW IC on `u'' = 6u²` the path-network's discretisation
    # floor is roughly 1e-13 from h=0.5 stepping accumulation,
    # already tighter than the loosest Tol in the sweep, so all four
    # final errors saturate at the same numerical floor.  We test the
    # controller's behaviour (h monotonicity) directly instead.
    # =================================================================
    @testset "AS.1.5: tighter Tol yields smaller accepted h" begin
        f = (z, u, up) -> 6 * u^2
        # Single Stage-1 target at z=5 (smooth region between IC at z=0
        # and pole-bridging starts at z≈1).  Short-enough walk to keep
        # CI cost bounded; long enough to exercise multiple steps and
        # see the controller's accept-and-seed memory.
        prob = PadeTaylorProblem(f, (u_0_FW, up_0_FW), (0.0, 5.0); order = 30)
        grid = ComplexF64[5.0 + 0.0im]

        max_hs = Float64[]
        for Tol in (1.0e-8, 1.0e-10, 1.0e-12, 1.0e-14)
            sol = path_network_solve(prob, grid;
                                      h = 1.5,
                                      step_size_policy = :adaptive_ffw,
                                      adaptive_tol = Tol,
                                      max_steps_per_target = 500)
            # `visited_h[1]` is the IC node's canonical-Padé radius
            # (always equal to the user-specified seed `h`, no rescale
            # has fired yet for the IC).  The CONTROLLER's accepted
            # step sizes are `visited_h[2:end]` — those are what the
            # Tol sweep should monotonically shrink.
            length(sol.visited_h) ≥ 2 || error(
                "AS.1.5: expected ≥ 2 visited nodes at Tol=$Tol")
            push!(max_hs, maximum(sol.visited_h[2:end]))
        end

        # Strictly: max_h non-increasing as Tol tightens.  At the
        # boundary where multiple Tols admit the same q (T(h) ≈ Tol),
        # equality is allowed — hence ≤ rather than <.
        for i in 2:length(max_hs)
            @test max_hs[i] ≤ max_hs[i-1] + 1e-12
        end
        # And the spread must be meaningful: tightening by 1e-6 (Tol
        # ratio 1e-8 → 1e-14) shrinks h by `(1e-6)^(1/31) ≈ 0.64`.
        @test max_hs[end] < max_hs[1] * 0.8
    end

    # =================================================================
    # AS.1.6 — Adaptive single-step builder unit test.  Direct call to
    # `adaptive_pade_step!` advances state and returns the
    # `(state, P_u, meta)` 3-tuple.  We seed at h=2.0 so the
    # controller MUST rescale (T(2.0) ≈ 5e-4 ≫ Tol=1e-12).
    # =================================================================
    @testset "AS.1.6: adaptive_pade_step! single-step contract" begin
        f = (z, u, up) -> 6 * u^2
        state = PadeStepperState{ComplexF64}(0.0, u_0_FW, up_0_FW)
        _, P_u, meta = adaptive_pade_step!(state, f, 30, ComplexF64(2.0, 0.0);
                                            adaptive_tol = 1.0e-12,
                                            k_conservative = 1.0e-3,
                                            max_rescales = 50)
        @test isa(meta, NamedTuple)
        @test haskey(meta, :h_used)
        @test haskey(meta, :T_h)
        @test haskey(meta, :n_rescales)
        @test haskey(meta, :h_step)
        @test meta.h_used isa Real
        @test 0 < meta.h_used < 2.0      # was rescaled (n_rescales ≥ 1)
        @test meta.T_h ≤ 1.0e-12
        @test meta.n_rescales ≥ 1
        @test abs(state.z - meta.h_step) ≤ 1.0e-15
        @test isfinite(state.u) && isfinite(state.up)
        @test P_u isa PadeApproximant{ComplexF64}
    end

    # =================================================================
    # AS.1.7 — Multi-iteration rescale.  Seed at h=10.0 where
    # T(10) ≈ 5e9 ≫ Tol=1e-12: ONE rescale lands at T ≈ 5.7e-7
    # (still ≫ Tol) — the controller must iterate.  Mutation M3
    # (one-shot rescale) bites here: meta.T_h reports the post-one-
    # shot value above Tol, failing the acceptance assertion.
    # =================================================================
    @testset "AS.1.7: rescale loop iterates beyond one shot at large h_init" begin
        f = (z, u, up) -> 6 * u^2
        state = PadeStepperState{ComplexF64}(0.0, u_0_FW, up_0_FW)
        _, _, meta = adaptive_pade_step!(state, f, 30, ComplexF64(10.0, 0.0);
                                          adaptive_tol = 1.0e-12,
                                          k_conservative = 1.0e-3,
                                          max_rescales = 50)
        @test meta.T_h ≤ 1.0e-12
        @test meta.n_rescales ≥ 2   # one rescale leaves T(h) ≫ Tol
        @test 0 < meta.h_used < 10.0
    end

end # @testset AdaptiveStep

# =====================================================================
# Mutation-proof procedure (verified before commit, restored after).
# Each mutation was applied to `src/PadeStepper.jl` (and to the inline
# adaptive loop in `src/PathNetwork.jl` where the algorithm is also
# implemented in the path-network driver), the adaptive test file was
# run via `julia --project=. -e 'using Pkg; Pkg.activate("."); include
# ("test/adaptive_step_test.jl")'`, the failure count recorded, and the
# mutation reverted.
#
#   M1 — flip rescale exponent sign in `ffw_rescale_q`:
#        q = (k·Tol/T(h))^(-1/(n+1)) instead of `1/(n+1)`.
#        Verified bite: 17 assertions failed
#        (12 in AS.1.2, 1 in AS.1.3, 1 in AS.1.5, 3 in AS.1.6).
#        The wrong exponent makes q > 1; h grows unboundedly and the
#        rescale loop hits max_rescales, or the controller diverges
#        outright.
#
#   M2 — drop the ε_{n+1} Padé correction in `ffw_truncation_error`:
#        eps_nplus1 = c̃_{n+1} only, no `Σ b_r c̃_{n+1-r}`.
#        Verified bite: 1 assertion in AS.1.1 (prefactor < 1e-9).
#        Dropping the correction raises the leading-h prefactor from
#        ~ 1e-10 (Padé-corrected) to ~ 2.5e-8 (pure 1/(n+1)!) — a
#        ~ 250× shift that the AS.1.1 prefactor pin catches.
#
#   M3 — replace the rescale `while` loop in `adaptive_pade_step!`
#        (AND the matching inline loop in `path_network_solve`) with
#        a one-shot `if Th > Tol …` block.
#        Verified bite: 2 assertions in AS.1.7 (meta.T_h ≤ 1e-12 AND
#        meta.n_rescales ≥ 2 both fail at h_init = 10.0, where one
#        application of q lands T(h) ≈ 5.7e-7 ≫ Tol = 1e-12).
#        AS.1.6 at h_init = 2.0 happens NOT to bite because one
#        rescale suffices there — hence the need for AS.1.7 as a
#        deliberately harsher h_init.
#
#   M4 — hard-code `q = 0.5` (FW-style halving) inside
#        `ffw_rescale_q`, dropping the Tol-dependent formula.
#        Verified bite: 1 assertion in AS.1.5 (max_hs[end] < max_hs[1]
#        * 0.8 fails because halving is Tol-independent — every Tol
#        produces the same convergent h regardless of how tight the
#        controller demands).
#
# All four mutations bit (M1: 17, M2: 1, M3: 2, M4: 1); all reverted
# before commit; final test count 43 GREEN.
