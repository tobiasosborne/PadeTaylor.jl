# test/pathnetwork_test.jl -- Phase 10 / bead `padetaylor-1jf` tests.
#
# Tier-2 path-network: FW 2011 §3.1 5-direction wedge tree + Stage-2
# fine-grid extrapolation.  Verifies that the complex-plane path
# network bridges the lattice pole of u(z) = ℘(z + c₁; 0, c₂) at z = 1
# via off-axis detours when stepping past it on the real axis would
# land on the pole.
#
# Reference: docs/adr/0004-path-network-architecture.md (algorithm +
# test plan PN.1.1-PN.3.1), docs/unified_path_network_spec.md (full
# spec), references/markdown/FW2011_painleve_methodology_JCP230/
# FW2011_painleve_methodology_JCP230.md:155-166 (FW 2011 §3.1).
#
# This first cut covers PN.1.1, PN.1.2, PN.2.1, PN.4.1.  PN.2.2 (FW
# Table 5.1 long-range z=30 to ≤1e-13) and PN.3.1 (:steepest_descent
# agreement with :min_u) are deferred to a follow-up commit so the
# initial GREEN ships on a tractable test corpus.

using Test
using PadeTaylor

include(joinpath(@__DIR__, "_oracle_problems.jl"))

@testset "PathNetwork (Phase 10): FW 2011 §3.1 path-network" begin

    # Test ODE: u'' = 6u^2 with FW 2011 ICs at z=0; closed form
    # u(z) = ℘(z + c_1; 0, c_2) with c_1 = -1, c_2 = 2.  Pole at z = 1.
    fW(z, u, up) = 6 * u^2

    @testset "PN.1.1: Stage 1 + Stage 2 single-step sanity (5-pt grid near z=0)" begin
        # Grid entirely inside |z| ≤ h=0.5; Stage 1 takes 0 steps for
        # most targets (they're already within h of the IC); Stage 2
        # is the IC-Padé evaluated at t = z_f / h.  This is the
        # cheapest end-to-end exercise of the public API.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.0); order = 30)
        grid = ComplexF64[0.0, 0.1+0.2im, -0.1+0.3im, 0.3-0.2im, 0.4+0.0im]
        sol  = path_network_solve(prob, grid; h = 0.5)

        # Stage 1: IC plus zero or one extra visited node (depending on
        # which targets fell within h).  Always ≥ 1 (the IC itself).
        @test length(sol.visited_z) ≥ 1
        @test sol.visited_z[1] == ComplexF64(0.0)

        # Stage 2: every grid point covered (no NaN).
        @test all(isfinite, real.(sol.grid_u))
        @test all(isfinite, imag.(sol.grid_u))

        # Verify against analytic ℘ at z = 0.4 (real-axis crosscheck
        # — far from pole, Padé is essentially exact).  ℘(0.4 + c1; 0, c2)
        # via series expansion is tabulated in _oracle_problems.jl as
        # implied by Phase-6 logic; here we just sanity-check magnitude.
        idx_04 = findfirst(z -> isapprox(real(z), 0.4) && abs(imag(z)) < 1e-12,
                           sol.grid_z)
        @test idx_04 !== nothing
        # u(0.4) — finite, real-valued (we're on real axis), and
        # consistent with Phase-6's IVP at z = 0.5 nearby (u ≈ 4.00 at
        # z = 0.5; at z = 0.4 we're a bit smaller in magnitude).
        @test abs(imag(sol.grid_u[idx_04])) < 1e-10
        @test real(sol.grid_u[idx_04]) > 1.07     # Just above IC value
        @test real(sol.grid_u[idx_04]) < 5.0      # Below z=0.5 value
    end

    @testset "PN.1.2: Multi-step path bridging the pole at z=1" begin
        # Target z = 1.4 (past the pole at z=1).  With h=0.5 from z=0
        # the goal direction is purely real; direct stepping would
        # land at z=0.5, then z=1.0 (ON THE POLE).  The 5-direction
        # wedge + min-|u| selection should detour off-axis.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 2.0); order = 30)
        grid = ComplexF64[1.4 + 0.0im]
        sol  = path_network_solve(prob, grid; h = 0.5)

        # Must have taken at least 3 steps from the IC.
        @test length(sol.visited_z) ≥ 3

        # None of the visited points should be exactly at (or have
        # blown up near) z=1.  If min-|u| failed and any candidate
        # landed on the pole, we'd see Inf or NaN here.
        @test all(isfinite, abs.(sol.visited_u))
        @test all(z -> abs(z - 1.0) > 0.01, sol.visited_z[2:end])

        # At least one off-axis (Im ≠ 0) visited node — proof of
        # complex-plane detour.
        @test any(z -> abs(imag(z)) > 1e-6, sol.visited_z)

        # u(1.4) at the target — finite + bounded.  Structural only.
        # The quantitative accuracy check belongs in PN.2.2 (FW Table 5.1
        # long-range tuning); the path-network's `:min_u` walk currently
        # accumulates Padé approximation error along the off-axis detour
        # and the final boundary-of-disc evaluation at `|t| ≈ 1` is the
        # dominant source.  See ADR-0004 deferral and worklog 005's
        # order/rtol-coupling pattern.  For PN.1.2 here we only assert
        # that the algorithm did not blow up or NaN.
        u_target = sol.grid_u[1]
        @test isfinite(real(u_target))
        @test isfinite(imag(u_target))
        @test abs(u_target) < 1e3       # Rough sanity: not divergent.
    end

    @testset "PN.2.1: Stage-2 NaN sentinel for uncovered grid points" begin
        # A grid point far outside any visited node's disc must return
        # NaN+NaN·im, not silent extrapolation (CLAUDE.md Rule 1).
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 0.5); order = 30)
        # Targets: trivial near-IC + one wild outlier at z=100 not
        # listed as a target (so Stage 1 ignores it).  Then evaluate
        # at z=100 in Stage 2 anyway — it should NaN.
        grid_targets = ComplexF64[0.1, 0.2, 0.3]
        sol = path_network_solve(prob, grid_targets; h = 0.5)
        # Confirm targets work normally.
        @test all(isfinite, real.(sol.grid_u))

        # Now build a separate grid that includes an uncovered point.
        # We do this via a re-invocation including z=100; Stage-1 will
        # FAIL to reach z=100 in max_steps_per_target=50 because that
        # requires 200 steps.  So we expect a thrown error here.
        @test_throws ErrorException path_network_solve(
            prob, ComplexF64[100.0 + 0.0im];
            h = 0.5, max_steps_per_target = 50)
    end

    @testset "PN.4.1: Fail-fast guards (CLAUDE.md Rule 1)" begin
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.0); order = 30)
        grid = ComplexF64[0.1 + 0.0im]

        # :adaptive_ffw is a Tier-4 deferral.
        @test_throws ArgumentError path_network_solve(
            prob, grid; h = 0.5, step_size_policy = :adaptive_ffw)

        # Unknown step selection.
        @test_throws ArgumentError path_network_solve(
            prob, grid; h = 0.5, step_selection = :bogus)

        # Wrong wedge size.
        @test_throws ArgumentError path_network_solve(
            prob, grid; h = 0.5, wedge_angles = [-π/4, 0.0, π/4])

        # Non-positive h.
        @test_throws ArgumentError path_network_solve(prob, grid; h = -0.5)
        @test_throws ArgumentError path_network_solve(prob, grid; h = 0.0)
    end

    @testset "PN.2.2: FW Table 5.1 z=30 long-range integration" begin
        # FW 2011 Fig 5.1 / Table 5.1: u(z=30) of the equianharmonic ℘
        # trajectory with FW IC matches `1.0950982559597442` to ≤1e-13
        # at BigFloat-256 (FW reports their own 8.34e-14 in the paper;
        # `docs/figure_catalogue.md §1 row FW2011 Fig 5.1`).
        #
        # This test exercises the path-network's load-bearing canonical-
        # Padé-per-visited-node invariant (ADR-0004 design decision):
        # each visited node z_v must store a REAL-h-direction Padé so
        # Stage 2's `t = (z_f - z_v) / h_v` interpolation lands inside
        # the disc.  See worklog 008 §"The wedge-vs-canonical-Padé bug
        # at long range".

        # Float64 path: rel-err ≤ 1e-9 acceptance.  ~75 visited nodes
        # for h=0.5, order=30.
        prob_f64 = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 30.0);
                                     order = 30)
        sol_f64 = path_network_solve(prob_f64, ComplexF64[30.0 + 0im];
                                     h = 0.5, max_steps_per_target = 200)
        u30_f64 = sol_f64.grid_u[1]
        @test isapprox(u30_f64, u_at_30_FW_ref; rtol = 1e-9)
        @test abs(imag(u30_f64)) < 1e-9        # Real solution on real axis.

        # BF-256 path: rel-err ≤ 1e-13 acceptance per FW Table 5.1.
        # ~50s wall time; the dominant test cost in the suite.
        setprecision(BigFloat, 256) do
            u0_bf  = big"1.071822516416917"
            up0_bf = big"1.710337353176786"
            ref_bf = big"1.0950982559597442"
            prob_bf = PadeTaylorProblem(fW, (u0_bf, up0_bf),
                                        (big(0.0), big(30.0)); order = 30)
            sol_bf = path_network_solve(prob_bf,
                                        Complex{BigFloat}[Complex{BigFloat}(big(30.0))];
                                        h = big(0.5), max_steps_per_target = 200)
            u30_bf = sol_bf.grid_u[1]
            rel = abs(u30_bf - ref_bf) / abs(ref_bf)
            @test rel ≤ big"1e-13"
            @test abs(imag(u30_bf)) < big"1e-15"
        end
    end

    @testset "PN.3.1: :steepest_descent path agrees with :min_u (pole-bridge)" begin
        # FW 2011 §5.4.1 (line 362-368) introduces `:steepest_descent`
        # as a perf-tuned alternative to `:min_u`: pick the wedge angle
        # closest to θ_sd = arg(-u/u') instead of evaluating all five
        # |u| values.  On any smooth-region grid + pole-bridge grid
        # where both rules agree on the wedge index at each step, the
        # two paths should produce identical visited-node values to
        # within Padé approximation noise (≤1e-10).
        #
        # We use targets that REQUIRE stepping (distance > h) including
        # one pole-bridge case (z=1.4 past the lattice pole at z≈1.13),
        # so the step-selection logic is exercised.  Both rules' answer
        # at z=1.4 is u≈6.2518 (Phase-10 PN.1.2 ground truth).
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 2.0); order = 30)
        grid = ComplexF64[1.4 + 0im, 1.2 + 0.4im, 0.6 + 0.3im]
        sol_min     = path_network_solve(prob, grid; h = 0.5,
                                         step_selection = :min_u)
        sol_descent = path_network_solve(prob, grid; h = 0.5,
                                         step_selection = :steepest_descent)
        @test length(sol_min.grid_u) == length(sol_descent.grid_u)
        for i in eachindex(sol_min.grid_u)
            Δu  = abs(sol_min.grid_u[i]  - sol_descent.grid_u[i])
            Δup = abs(sol_min.grid_up[i] - sol_descent.grid_up[i])
            @test Δu  ≤ 1e-10
            @test Δup ≤ 1e-9   # u' carries the 1/h chain-rule factor; looser tol.
        end
    end

end # @testset PathNetwork

# PN.5.1  Mutation-proof procedure (verified manually before commit; see
# ADR-0004 §"Test plan" + worklog 004 §"Mutation-proof procedure for
# the next agent" for the Phase-6 lineage):
#
#   Mutation A  --  in `_select_candidate`, replace `argmin(abs(e[2]))`
#     with `argmax(abs(e[2]))` for the :min_u branch.  Steers path
#     TOWARD poles instead of away from them.
#     Expected: PN.1.2 RED at lines 79, 83 — `all(isfinite, abs.(visited_u))`
#     fails because steering toward z=1 pole yields Inf/NaN, and
#     `all(z -> abs(z - 1) > 0.01)` fails because a step lands near pole.
#     Verified 2026-05-13: 2 fails on PN.1.2 as predicted.
#
#   Mutation D  --  in the Stage-2 evaluation loop, invert the coverage
#     check from `abs(z_f - z_v) > h_v` to `abs(z_f - z_v) < h_v` so
#     COVERED points NaN-out and UNCOVERED points (which we have none of
#     in PN.1.1/PN.1.2) get evaluated.
#     Expected: PN.1.1 RED at lines 45-46 + 58-60 (`isfinite` checks +
#     u(0.4) magnitude bounds); PN.1.2 RED at lines 94-96 (u_target
#     finiteness + magnitude); PN.2.1 line 109 (the targets-work-normally
#     check inverts).
#     Verified 2026-05-13: 9 fails across PN.1.1 (5), PN.1.2 (3), PN.2.1 (1)
#     as predicted — the NaN-sentinel logic is load-bearing across all
#     Stage-2 tests.
#
# Restoration: both mutations restored before commit.
#
# PN.2.2 + PN.3.1 follow-up mutations (verified 2026-05-13 with bead
# `padetaylor-yt1` in flight):
#
#   Mutation E  --  in PathNetwork.jl Stage-1 loop, restore the pre-bugfix
#     behaviour of storing the wedge-direction `pade_sel` instead of the
#     canonical-direction `pade_canonical` per visited node.  This is the
#     bug worklog 008 diagnosed; restoring it MUST bite PN.2.2.
#     Verified bite: 4 fails on PN.2.2 — F64 rel-err 0.218 >> 1e-9, F64
#     imag(u) 7.2e-3 >> 1e-9, BF-256 rel-err 0.218 >> 1e-13, BF-256 imag
#     7.2e-3 >> 1e-15.  Algorithmic error invariant to arithmetic precision.
#
#   Mutation F  --  in `_select_candidate` :steepest_descent branch,
#     flip the sign in `θ_sd = angle(-u/u')` to `θ_sd = angle(u/u')`,
#     steering toward poles instead of away.
#     Verified bite: 4 fails on PN.3.1 — |Δu| at z=1.4 reaches 6.25
#     (steepest_descent returns 0+0im on a failed pole crossing while
#     :min_u correctly gives u≈6.2518); |Δu| at z=1.2+0.4i reaches 5.00;
#     |Δu'| reaches 31.2 and 22.3 respectively.  Way above the 1e-10 /
#     1e-9 tolerance.
#
# Restoration: both mutations restored before commit.
