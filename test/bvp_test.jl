# test/bvp_test.jl -- Phase 11 / bead `padetaylor-804` tests.
#
# Tier-3 Chebyshev-Newton spectral-collocation BVP solver.  Validates
# against the DMSUITE-based Octave oracle at external/probes/bvp-oracle/
# (Octave 8.4.0, captured 2026-05-13).  Five groups mirror the oracle's
# five groups:
#
#   1. D₂ matrix sanity (chebdif primitives)
#   2. Linear BVP u''=u on [-1,1] (closed form cosh(t)/cosh(1))
#   3. Nonlinear PI BVP on real segment [-18,-14] (FW 2011 Fig 4.6 case)
#   4. Nonlinear PI BVP on real segment [0.1,0.5] (mpmath-pinned BCs)
#   5. Barycentric interpolation at intermediate points
#
# Plus: callable returns (u, u'); fail-fast guards (BV.4.x); BigFloat-256
# sanity (BV.5.1).
#
# Reference: references/bvp_recipe.md (canonical algorithm),
# docs/adr/0004-path-network-architecture.md (Tier-3 plan),
# references/markdown/FW2011_painleve_methodology_JCP230/
# FW2011_painleve_methodology_JCP230.md:176-200 (FW 2011 §3.2 spec).

using Test
using PadeTaylor
using PadeTaylor: BVP   # for access to private helpers in BV.1.1

include(joinpath(@__DIR__, "_oracle_bvp.jl"))

@testset "BVP (Phase 11): Chebyshev-Newton spectral collocation" begin

    # PI ODE: u'' = 6u^2 + z.  Analytic Jacobian: ∂F/∂u = 12u.
    f_PI(z, u)    = 6 * u^2 + z
    ∂f_PI(z, u)   = 12 * u

    # Linear ODE u'' = u for sanity.
    f_lin(z, u)   = u
    ∂f_lin(z, u)  = one(u)

    @testset "BV.1.1: D₂ matrix entries match DMSUITE (chebdif primitives)" begin
        # N=4: full 5×5 D₁ and D₂.  Tolerance 1e-12 — both sides Float64;
        # the only floating arithmetic difference is operation ordering.
        nodes_N4_jl = [cos(j * π / 4) for j in 0:4]
        @test isapprox(nodes_N4_jl, test_d2_N4_nodes; atol = 1e-14)

        D1_N4 = BVP._chebyshev_D1(nodes_N4_jl, Float64, 4)
        @test isapprox(D1_N4, test_d2_N4_D1; atol = 1e-12)
        @test isapprox(D1_N4 * D1_N4, test_d2_N4_D2; atol = 1e-11)

        # N=8: full 9×9.
        nodes_N8_jl = [cos(j * π / 8) for j in 0:8]
        D1_N8 = BVP._chebyshev_D1(nodes_N8_jl, Float64, 8)
        @test isapprox(D1_N8, test_d2_N8_D1; atol = 1e-12)
        @test isapprox(D1_N8 * D1_N8, test_d2_N8_D2; atol = 1e-10)

        # N=16: partial — only diagonal + super-diagonal of D₂ pinned in oracle.
        nodes_N16_jl = [cos(j * π / 16) for j in 0:16]
        D1_N16 = BVP._chebyshev_D1(nodes_N16_jl, Float64, 16)
        D2_N16 = D1_N16 * D1_N16
        @test isapprox([D2_N16[i,i] for i in 1:17], test_d2_N16_D2_diag; atol = 1e-8)
        @test isapprox([D2_N16[i, i+1] for i in 1:16], test_d2_N16_D2_supdiag; atol = 1e-8)
    end

    @testset "BV.1.2: Linear BVP u'' = u matches cosh(t)/cosh(1)" begin
        # N=8: 3 Newton iters per oracle.  Solution error ~7e-10 is the
        # Chebyshev spectral truncation floor, not Newton failure.
        sol8 = bvp_solve(f_lin, ∂f_lin, -1.0, 1.0, 1.0, 1.0; N = 8)
        @test sol8.iterations ≤ 3
        @test sol8.residual_inf ≤ 1e-13
        @test isapprox(sol8.u_nodes, test_bvp_lin_N8_u; atol = 1e-12)
        # Versus the closed-form cosh(t)/cosh(1) — should saturate at
        # the Chebyshev N=8 truncation floor (~1e-9).
        ref = [cosh(t) / cosh(1.0) for t in sol8.nodes_t]
        @test maximum(abs, sol8.u_nodes .- ref) < 1e-8

        # N=16: 1 Newton iter (linear problem converges in one step).
        sol16 = bvp_solve(f_lin, ∂f_lin, -1.0, 1.0, 1.0, 1.0; N = 16)
        @test sol16.iterations ≤ 2
        @test sol16.residual_inf ≤ 1e-12
        @test isapprox(sol16.u_nodes, test_bvp_lin_N16_u; atol = 1e-11)
        ref16 = [cosh(t) / cosh(1.0) for t in sol16.nodes_t]
        @test maximum(abs, sol16.u_nodes .- ref16) < 1e-13
    end

    @testset "BV.1.3: Nonlinear PI BVP on real segment [-18,-14]" begin
        # FW 2011 Fig 4.6 case.  IC: u ≈ √(-z/6) asymptotic.
        sqrt_guess(z) = sqrt(-z / 6)
        sol = bvp_solve(f_PI, ∂f_PI, -18.0, -14.0,
                        sqrt_guess(-18.0), sqrt_guess(-14.0);
                        N = 20, initial_guess = sqrt_guess)
        @test sol.iterations ≤ 4
        @test sol.residual_inf ≤ 1e-11
        # Solution values at all 21 collocation nodes match oracle.
        @test isapprox(sol.u_nodes, test_bvp_pi_real_N20_u; atol = 1e-10)
        # Affine-mapped z nodes match oracle.
        @test isapprox(sol.nodes_z, test_bvp_pi_real_N20_nodes_z; atol = 1e-12)
    end

    @testset "BV.1.4: Nonlinear PI BVP on real segment [0.1,0.5]" begin
        # BCs are mpmath-pinned values from the Octave oracle's PI solve.
        # CRITICAL: u'' = 6u^2 + z (PI), NOT u'' = 6u^2 (equianharmonic ℘).
        # See worklog-note in test/_oracle_bvp.jl Group 4 header.
        sol = bvp_solve(f_PI, ∂f_PI, 0.1, 0.5,
                        test_bvp_pi_complex_N24_bc_a,
                        test_bvp_pi_complex_N24_bc_b;
                        N = 24)
        @test sol.iterations ≤ 5
        @test sol.residual_inf ≤ 1e-11
        @test isapprox(sol.u_nodes, test_bvp_pi_complex_N24_u; atol = 1e-10)
    end

    @testset "BV.2.1: Barycentric evaluation matches chebint oracle" begin
        # Linear N=8: eval at t* ∈ {-0.7, -0.3, 0, 0.4, 0.8}.
        # Since our affine map sends t* directly to z* = t* on [-1, 1]:
        sol8 = bvp_solve(f_lin, ∂f_lin, -1.0, 1.0, 1.0, 1.0; N = 8)
        for (i, t_star) in enumerate(test_baryeval_lin_N8_t_star)
            (u_at, _) = sol8(t_star)        # z* = t* on [-1, 1]
            @test isapprox(u_at, test_baryeval_lin_N8_p[i]; atol = 1e-12)
        end

        # PI real N=20: eval at z* values from oracle (under affine map).
        sqrt_guess(z) = sqrt(-z / 6)
        sol20 = bvp_solve(f_PI, ∂f_PI, -18.0, -14.0,
                          sqrt_guess(-18.0), sqrt_guess(-14.0);
                          N = 20, initial_guess = sqrt_guess)
        for (i, z_star) in enumerate(test_baryeval_pi_real_N20_z_star)
            (u_at, _) = sol20(z_star)
            @test isapprox(u_at, test_baryeval_pi_real_N20_p[i]; atol = 1e-10)
        end
    end

    @testset "BV.3.1: Callable returns (u, u') with correct endpoint derivatives" begin
        sqrt_guess(z) = sqrt(-z / 6)
        sol = bvp_solve(f_PI, ∂f_PI, -18.0, -14.0,
                        sqrt_guess(-18.0), sqrt_guess(-14.0);
                        N = 20, initial_guess = sqrt_guess)

        # At z = z_b = -14 (t = +1) — should recover BC and oracle derivative.
        (u_at_zb, up_at_zb) = sol(-14.0)
        @test isapprox(u_at_zb, sqrt_guess(-14.0); atol = 1e-10)
        @test isapprox(up_at_zb, test_bvp_pi_real_N20_dup_at_zb; atol = 1e-9)

        # At z = z_a = -18 (t = -1) — same.
        (u_at_za, up_at_za) = sol(-18.0)
        @test isapprox(u_at_za, sqrt_guess(-18.0); atol = 1e-10)
        @test isapprox(up_at_za, test_bvp_pi_real_N20_dup_at_za; atol = 1e-9)
    end

    @testset "BV.4.1: Fail-fast guards (CLAUDE.md Rule 1)" begin
        # N too small.
        @test_throws ArgumentError bvp_solve(f_lin, ∂f_lin, -1.0, 1.0, 1.0, 1.0; N = 3)

        # maxiter zero.
        @test_throws ArgumentError bvp_solve(f_lin, ∂f_lin, -1.0, 1.0, 1.0, 1.0;
                                             N = 8, maxiter = 0)

        # Negative tol.
        @test_throws ArgumentError bvp_solve(f_lin, ∂f_lin, -1.0, 1.0, 1.0, 1.0;
                                             N = 8, tol = -1.0)

        # Non-convergent Newton (deliberately tight tol + wild BC).
        # The linear u''=u with BC u(±1) = 1e10 should require many iters but
        # since it's linear it converges in 1 step.  Make it non-convergent by
        # using a wildly-off initial_guess + maxiter = 1 on PI (nonlinear).
        bad_init(z) = 1000.0 + z   # very far from any PI solution
        @test_throws ErrorException bvp_solve(f_PI, ∂f_PI, -18.0, -14.0,
                                              sqrt(18 / 6), sqrt(14 / 6);
                                              N = 8, maxiter = 1,
                                              initial_guess = bad_init)

        # Out-of-segment evaluation throws DomainError.
        sol = bvp_solve(f_lin, ∂f_lin, -1.0, 1.0, 1.0, 1.0; N = 8)
        @test_throws DomainError sol(5.0)
    end

    @testset "BV.5.1: BigFloat-256 precision sanity (linear BVP)" begin
        # Run u''=u on [-1, 1] at BF-256.  At N=8 the Chebyshev truncation
        # floor (~1e-9) saturates BEFORE BF-256 helps — worklog-005 order/
        # precision-coupling pattern.  At N=16 BF-256 starts to overtake.
        setprecision(BigFloat, 256) do
            sol = bvp_solve(f_lin, ∂f_lin, big(-1.0), big(1.0), big(1.0), big(1.0); N = 16)
            @test sol.residual_inf ≤ big(1e-60)   # BF-256 floor much tighter than Float64
            ref = [cosh(t) / cosh(big(1.0)) for t in sol.nodes_t]
            # At N=16, Chebyshev truncation error for cosh on [-1,1] is ~1e-16.
            # BF-256 fully saturates that floor.
            @test maximum(abs, sol.u_nodes .- ref) < big(1e-15)
        end
    end

end # @testset BVP

# BV.6.1  Mutation-proof procedure (verified 2026-05-13 before commit;
# see references/bvp_recipe.md §4 for the Jacobian factor's derivation):
#
#   Mutation A  --  in `bvp_solve`, drop the (1/4) affine-scale factor,
#     using `scale = (z_b - z_a)^2` (4× too large).  Both residual and
#     Jacobian rescale incorrectly.
#     Verified bite: 19 fails across BV.1.2 (8/8), BV.1.3 (1/4), BV.1.4
#     (2/3), BV.2.1 (8/9), BV.3.1 (2/4), BV.5.1 (1/2).  The (z_b-z_a)²/4
#     factor is THE single most load-bearing element of the recipe — its
#     absence corrupts every nontrivial test.  This is the spec-drift
#     pattern the FW2011 line 190 "the Jacobian is trivial to calculate"
#     phrasing invites; the worklog-005 / Octave-oracle subagent already
#     caught one variant of it.
#
#   Mutation B  --  in `bvp_solve`, swap `u_a` and `u_b` in BOTH the
#     direct BC enforcement (`u_nodes[1] = u_a` instead of `u_b`) and
#     the precomputed `bc_col` contribution from the BC columns.
#     Linear BVP with symmetric BCs (u(±1) = 1) doesn't bite, but every
#     asymmetric test does.
#     Verified bite: 10 fails across BV.1.3 (1/4), BV.1.4 (2/3), BV.2.1
#     (4/9), BV.3.1 (3/4).
#
#   Mutation C  --  in `_barycentric_eval`, drop the `(-1)^j` alternating
#     sign (`sign_w = T(1)` always).  Newton solve unaffected (D₁/D₂
#     don't use barycentric weights), but the callable's interpolation
#     is wrong.
#     Verified bite: 8 fails, all in BV.2.1's eval-at-intermediate-points.
#
# Restoration: all three mutations restored before commit per CLAUDE.md
# Rule 4.  Matches the Phase-6 + Phase-10 inline mutation-proof pattern.
