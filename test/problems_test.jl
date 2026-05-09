# test/problems_test.jl -- Phase 6 tests.
#
# The Padé-vs-Taylor pole-bridge demonstration: ONE Padé built at z = 0
# with h_max = 1.5 spans the lattice pole of u(z) = ℘(z + c₁; 0, c₂)
# at z = 1.  In the rescaled variable t = z / h_max the pole sits at
# t = 1/1.5 ≈ 0.667, strictly INSIDE the segment [0, 1].  The same Padé,
# evaluated at t = z/1.5 for z ∈ [0, 1.5], must therefore give correct
# values both BEFORE the pole (z = 0.5, 0.95) and AFTER (z = 1.05, 1.4)
# — bridging the pole is the central claim of the analytic-continuation
# advantage of Padé over plain Taylor.
#
# Reference: docs/worklog/004-phase-6-pivot.md (full spec, oracle plan,
# v1/v2 acceptance scope).  Oracles cross-checked in
# external/probes/problems-oracle/verify.jl.
#
# v1 acceptance: tests 6.1.1–6.1.6 GREEN at the listed tolerances.  v2
# (FW Table 5.1 long-range integration to z = 30) is deferred — see
# bead padetaylor-8cr.

using Test
using PadeTaylor

include(joinpath(@__DIR__, "_oracle_problems.jl"))

@testset "Problems (Phase 6): Padé bridges the lattice pole at z=1" begin
    fW(z, u, up) = 6 * u^2

    @testset "6.1.1  sol(0.5) ~ WeierstrassP(0.5) to 1e-13 [3-source]" begin
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
        sol = solve_pade(prob; h_max = 1.5)
        u, up = sol(0.5)
        @test isapprox(u,  u_at_0_5;  rtol = 1e-13)
        @test isapprox(up, up_at_0_5; rtol = 1e-13)
    end

    @testset "6.1.2  sol(0.95) ~ WeierstrassP(0.95) to 1e-9 [near-pole, 3-source]" begin
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
        sol = solve_pade(prob; h_max = 1.5)
        u, up = sol(0.95)
        @test isapprox(u,  u_at_0_95;  rtol = 1e-9)
        @test isapprox(up, up_at_0_95; rtol = 1e-9)
    end

    @testset "6.1.3  sol(1.05) ~ WeierstrassP(1.05) to 1e-9 -- Padé bridges pole" begin
        # Headline: same Padé, evaluated at t = 1.05/1.5 = 0.7, i.e. just
        # PAST the pole at t = 0.667.  Plain Taylor diverges here (test 6.1.5);
        # Padé tracks the closed-form to 1e-9 because the denominator zero
        # at t ≈ 0.667 captures the lattice pole structure.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
        sol = solve_pade(prob; h_max = 1.5)
        u, up = sol(1.05)
        @test isapprox(u,  u_at_1_05;  rtol = 1e-9)
        @test isapprox(up, up_at_1_05; rtol = 1e-9)
    end

    @testset "6.1.4  sol(1.4) ~ WeierstrassP(1.4) to 1e-7 -- further past pole" begin
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
        sol = solve_pade(prob; h_max = 1.5)
        u, up = sol(1.4)
        @test isapprox(u,  u_at_1_4;  rtol = 1e-7)
        @test isapprox(up, up_at_1_4; rtol = 1e-7)
    end

    @testset "6.1.5  Padé wins, plain Taylor diverges -- the headline demo" begin
        # Same Taylor coefficients on both sides; the only difference is
        # whether we Padé-convert or truncate.  Past the natural radius of
        # convergence at z = 1, plain Taylor explodes by orders of magnitude;
        # Padé continues to track the closed-form ℘.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
        sol = solve_pade(prob; h_max = 1.5)
        u_pade, _ = sol(1.05)

        coefs_u = taylor_coefficients_2nd(fW, 0.0, u_0_FW, up_0_FW, 30)
        u_taylor = PadeTaylor.taylor_eval(coefs_u, 1.05)

        # Padé wins: rel err < 1e-9.
        @test isapprox(u_pade, u_at_1_05; rtol = 1e-9)
        # Taylor loses: relative error > 0.1 (truncation diverges past
        # natural radius of convergence at z = 1).
        @test abs(u_taylor - u_at_1_05) / abs(u_at_1_05) > 0.1
    end

    @testset "6.1.6  BigFloat-256: sol(1.05) at high precision (order = 40)" begin
        # Why order = 40, not 30: empirical sweep at z = 1.05 (just past the
        # pole) shows the (15, 15) diagonal Padé from order = 30 has an
        # inherent rational-approximation-error floor of ~3.5e-10 at
        # t = 0.7 — at THAT order BF-256 adds no value over Float64 because
        # both saturate at the same Padé-error floor.  Bumping to order = 40
        # changes the picture: Float64 *degrades* to ~1.6e-8 (SVD ill-
        # conditioning of the 41-by-41 Toeplitz block), while BF-256
        # cleanly converges to ~5e-15 (the IC-precision floor for the
        # 16-digit FW ICs after amplification through the pole).  This
        # makes 6.1.6 a meaningful test of the BF-256 path: at this
        # order, BF-256 is REQUIRED, not optional.
        setprecision(BigFloat, 256) do
            u0_bf  = BigFloat("1.071822516416917")
            up0_bf = BigFloat("1.710337353176786")
            zstart = BigFloat("0")
            zend   = BigFloat(3) / BigFloat(2)   # exact 1.5

            prob = PadeTaylorProblem(fW, (u0_bf, up0_bf), (zstart, zend); order = 40)
            sol = solve_pade(prob; h_max = BigFloat(3) / BigFloat(2))

            u_ref  = parse(BigFloat, u_at_1_05_80dps_str)
            up_ref = parse(BigFloat, up_at_1_05_80dps_str)

            u_eval, up_eval = sol(BigFloat(105) / BigFloat(100))
            @test isapprox(u_eval,  u_ref;  rtol = BigFloat("1e-13"))
            @test isapprox(up_eval, up_ref; rtol = BigFloat("1e-13"))
        end
    end
end

# 6.1.7  Mutation-proof procedure (verified manually before commit; see
# worklog 004 §"Mutation-proof procedure for the next agent"):
#
#   Mutation A  --  in `pade_step_with_pade!` (or `solve_pade`), replace
#     the `robust_pade(c̃, m, n)` call with one that returns the Taylor-
#     truncation Padé (numerator = c̃, denominator = [1.0]).
#     Expected: 6.1.3 + 6.1.4 RED (Taylor diverges past the pole at t=0.667).
#
#   Mutation B  --  in `(sol::PadeTaylorSolution)(z)`, return the
#     segment-start state `sol.y[k]` instead of evaluating the Padé at
#     the rescaled t.
#     Expected: 6.1.1 RED  (z = 0.5 returns u(0) ≈ 1.072 not u(0.5) ≈ 4.004).
