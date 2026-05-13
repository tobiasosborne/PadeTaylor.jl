# test/dispatcher_test.jl -- Phase 12 / bead `padetaylor-8lk` tests.
#
# Tier-3 Dispatcher: 1D ordered-chain composition of PathNetwork (IVP)
# and BVP (Chebyshev-Newton) segments per FW 2011 §4.4
# (`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:249-261`).
# Junction derivative-match diagnostic per the same paper, line 192
# ("typically set to 10⁻⁷ or 10⁻⁸").
#
# v1 scope: linear ordered chain of IVPSegment / BVPSegment specs.
#   - DP.1.1: single :ivp segment ≡ direct path_network_solve.
#   - DP.1.2: single :bvp segment ≡ direct bvp_solve.
#   - DP.2.1: two-segment IVP→BVP with junction derivative-match.
#   - DP.3.1: three-segment IVP→BVP→IVP (FW 2011 §4.4 composition).
#   - DP.4.1: fail-fast guards.
#   - DP.5.1: mutation-proof (procedure documented at bottom).
#
# v2 (deferred bead `padetaylor-k31`): 2D lattice dispatcher with
# automatic 5-point Laplacian edge detection per FW2011...md:204-208.
#
# Test ODE family: linear u'' = u (closed form u(t) = A·cosh(t) + B·sinh(t))
# is the spine.  IC u(0)=1, u'(0)=0 gives u(t)=cosh(t), u'(t)=sinh(t).
# Clean closed-form cross-checks for both IVP and BVP branches, non-trivial
# smooth solution for junction tests.  The Phase-10 equianharmonic ℘
# setup (u'' = 6u^2 with FW 2011 ICs at z=0) is used in DP.1.1 to
# verify the passthrough adapter on a problem the existing path-network
# tests already pin.

using Test
using PadeTaylor
using PadeTaylor: Dispatcher

include(joinpath(@__DIR__, "_oracle_problems.jl"))

@testset "Dispatcher (Phase 12): FW 2011 §4.4 IVP↔BVP composition" begin

    # Linear ODE u'' = u (the spine; closed form u(t) = cosh(t) at IC
    # u(0)=1, u'(0)=0).
    f_lin_2(z, u, up) = u            # PathNetwork: 2nd-order f(z, u, u')
    f_lin_1(z, u)     = u            # BVP: 1st-order f(z, u)
    ∂f_lin_1(z, u)    = one(u)

    # Equianharmonic ℘ ODE u'' = 6u^2 (Phase-10 PathNetwork test bed).
    fW_2(z, u, up) = 6 * u^2

    @testset "DP.1.1: Single IVPSegment passthrough ≡ path_network_solve" begin
        # Dispatcher with a lone :ivp segment must equal a direct
        # path_network_solve call on the same problem + grid.  Verifies
        # the IVP-side adapter is pure passthrough.
        prob = PadeTaylorProblem(fW_2, (u_0_FW, up_0_FW), (0.0, 1.0); order = 30)
        dense_grid = ComplexF64[0.1, 0.2, 0.3, 0.4]
        seg = IVPSegment(ComplexF64(0.5), dense_grid)

        sol = dispatch_solve(prob, fW_2, (z, u) -> 12 * u, [seg]; h = 0.5)
        sol_direct = path_network_solve(prob, dense_grid; h = 0.5)

        @test sol.grid_z ≈ sol_direct.grid_z
        @test sol.grid_u ≈ sol_direct.grid_u
        @test sol.grid_up ≈ sol_direct.grid_up
        @test all(sol.grid_region .== :ivp)
        @test length(sol.bvp_solutions) == 0
        @test length(sol.junction_z) == 0
        @test length(sol.ivp_solutions) == 1
        @test sol.segment_kinds == [:ivp]
    end

    @testset "DP.1.2: Single BVPSegment passthrough ≡ bvp_solve" begin
        # Dispatcher with a lone :bvp segment must equal a direct
        # bvp_solve call.  Linear ODE u''=u on [-1, 1] with u(±1) = 1
        # — the Phase-11 BV.1.2 reference setup.  prob.y0[1] = 1.0
        # is taken as the BVP's left BC (u_a).
        prob = PadeTaylorProblem(f_lin_2, (1.0, 0.0), (-1.0, 1.0); order = 30)
        dense_grid = ComplexF64[-0.5, 0.0, 0.5]
        seg = BVPSegment(ComplexF64(1.0), ComplexF64(1.0), dense_grid)

        sol = dispatch_solve(prob, f_lin_1, ∂f_lin_1, [seg];
                             h = 0.5, N_bvp = 8)
        sol_direct = bvp_solve(f_lin_1, ∂f_lin_1, -1.0, 1.0, 1.0, 1.0; N = 8)

        @test sol.bvp_solutions[1].iterations == sol_direct.iterations
        @test sol.bvp_solutions[1].u_nodes ≈ sol_direct.u_nodes
        @test sol.segment_kinds == [:bvp]
        @test all(sol.grid_region .== :bvp)
        @test length(sol.junction_z) == 0
        # Dense-grid eval matches barycentric eval of the direct sol.
        for (i, z) in enumerate(dense_grid)
            (u_ref, _) = sol_direct(z)
            @test isapprox(sol.grid_u[i], u_ref; atol = 1e-12)
        end
    end

    @testset "DP.2.1: IVP→BVP junction with derivative-match diagnostic (FW 2011 line 192)" begin
        # Linear ODE u'' = u, IC u(0)=1, u'(0)=0 → u(t)=cosh(t), u'(t)=sinh(t).
        # IVP from 0 to 0.3, then BVP from 0.3 to 0.6 with u_a from
        # IVP terminus and u_b = cosh(0.6).  Junction at z=0.3.
        # The BVP recovers u'_BVP(0.3) via barycentric on D₁·u_nodes;
        # FW2011 line 192 expects |u'_IVP - u'_BVP| ≤ 1e-7 typical.
        prob = PadeTaylorProblem(f_lin_2, (1.0, 0.0), (0.0, 0.3); order = 30)
        ivp_grid = ComplexF64[0.1, 0.2]
        bvp_grid = ComplexF64[0.35, 0.4, 0.45, 0.5, 0.55]
        u_b_target = ComplexF64(cosh(0.6))

        seg_ivp = IVPSegment(ComplexF64(0.3), ivp_grid)
        seg_bvp = BVPSegment(ComplexF64(0.6), u_b_target, bvp_grid)

        sol = dispatch_solve(prob, f_lin_1, ∂f_lin_1,
                             [seg_ivp, seg_bvp];
                             h = 0.5, N_bvp = 16,
                             derivative_match_tol = 1e-7)

        # One interior junction at z=0.3.
        @test length(sol.junction_z) == 1
        @test sol.junction_z[1] ≈ ComplexF64(0.3)
        # Δu is zero by construction (BC enforced from IVP terminus).
        # Δu' is the diagnostic; linear ODE + N=16 should be ≤ 1e-9.
        (Δu, Δup) = sol.junction_match[1]
        @test Δu ≤ 1e-12
        @test Δup ≤ 1e-7

        # Verify against closed form at sample points within each region.
        # cosh at z=0.2 (IVP region):
        idx_02 = findfirst(z -> isapprox(real(z), 0.2), sol.grid_z)
        @test idx_02 !== nothing
        @test isapprox(sol.grid_u[idx_02], cosh(0.2); atol = 1e-12)
        @test sol.grid_region[idx_02] == :ivp
        # cosh at z=0.45 (BVP region):
        idx_045 = findfirst(z -> isapprox(real(z), 0.45), sol.grid_z)
        @test idx_045 !== nothing
        @test isapprox(sol.grid_u[idx_045], cosh(0.45); atol = 1e-10)
        @test sol.grid_region[idx_045] == :bvp
    end

    @testset "DP.3.1: Three-segment IVP→BVP→IVP (FW 2011 §4.4 composition pattern)" begin
        # IVP from 0 to 0.3, BVP from 0.3 to 0.6, IVP from 0.6 to 0.9.
        # Both junctions get derivative-match diagnostics; final state
        # at z≈0.75 in the third IVP segment should match cosh(0.75)
        # within the cumulative IVP+BVP+IVP error budget.
        prob = PadeTaylorProblem(f_lin_2, (1.0, 0.0), (0.0, 0.3); order = 30)
        seg1 = IVPSegment(ComplexF64(0.3), ComplexF64[0.15])
        seg2 = BVPSegment(ComplexF64(0.6),
                          ComplexF64(cosh(0.6)),
                          ComplexF64[0.45])
        seg3 = IVPSegment(ComplexF64(0.9), ComplexF64[0.75])

        sol = dispatch_solve(prob, f_lin_1, ∂f_lin_1,
                             [seg1, seg2, seg3];
                             h = 0.4, N_bvp = 16,
                             derivative_match_tol = 1e-7)

        # Two interior junctions.
        @test length(sol.junction_z) == 2
        @test sol.junction_z[1] ≈ ComplexF64(0.3)
        @test sol.junction_z[2] ≈ ComplexF64(0.6)
        @test sol.junction_match[1][1] ≤ 1e-12   # |Δu| at first junction
        @test sol.junction_match[1][2] ≤ 1e-7    # |Δu'|
        @test sol.junction_match[2][1] ≤ 1e-12   # |Δu| at second junction
        @test sol.junction_match[2][2] ≤ 1e-7    # |Δu'|

        # Sanity: third IVP segment's grid_u contains valid cosh values.
        # Cumulative error budget: IVP (small, ~1e-10) + BVP (1e-10) +
        # IVP from BVP-derived IC (1e-7 inherited).  Loose 1e-5 tol.
        idx_075 = findfirst(z -> isapprox(real(z), 0.75), sol.grid_z)
        @test idx_075 !== nothing
        @test isapprox(sol.grid_u[idx_075], cosh(0.75); atol = 1e-5)
        @test sol.grid_region[idx_075] == :ivp

        @test sol.segment_kinds == [:ivp, :bvp, :ivp]
        @test length(sol.ivp_solutions) == 2
        @test length(sol.bvp_solutions) == 1
    end

    @testset "DP.4.1: Fail-fast guards (CLAUDE.md Rule 1)" begin
        prob = PadeTaylorProblem(f_lin_2, (1.0, 0.0), (0.0, 1.0); order = 30)
        seg_ivp = IVPSegment(ComplexF64(0.5), ComplexF64[0.2])
        seg_bvp = BVPSegment(ComplexF64(0.5), ComplexF64(1.0), ComplexF64[0.3])

        # Empty segments.
        @test_throws ArgumentError dispatch_solve(
            prob, f_lin_1, ∂f_lin_1, IVPSegment{Float64}[]; h = 0.5)

        # Negative derivative_match_tol.
        @test_throws ArgumentError dispatch_solve(
            prob, f_lin_1, ∂f_lin_1, [seg_ivp]; h = 0.5,
            derivative_match_tol = -1.0)

        # Negative h.
        @test_throws ArgumentError dispatch_solve(
            prob, f_lin_1, ∂f_lin_1, [seg_ivp]; h = -0.5)

        # N_bvp < 4 (delegated to bvp_solve; we either pre-validate or
        # let the underlying call throw — either is acceptable per Rule 1).
        @test_throws ArgumentError dispatch_solve(
            prob, f_lin_1, ∂f_lin_1, [seg_bvp]; h = 0.5, N_bvp = 3)

        # Strict mode with a deliberately-tight tol smaller than the
        # spectral floor (~ N² · eps), then a tiny N_bvp.  At N_bvp = 4
        # on a [0, 0.6] segment with u_a, u_b = cosh-values, the BVP's
        # barycentric u' will deviate from the IVP's analytic u' by
        # ≳ 1e-6, exceeding tol = 1e-12.  Strict-mode → ErrorException.
        seg_ivp_tight = IVPSegment(ComplexF64(0.3), ComplexF64[0.15])
        seg_bvp_tight = BVPSegment(ComplexF64(0.6),
                                   ComplexF64(cosh(0.6)),
                                   ComplexF64[0.45])
        @test_throws ErrorException dispatch_solve(
            prob, f_lin_1, ∂f_lin_1, [seg_ivp_tight, seg_bvp_tight];
            h = 0.5, N_bvp = 4,
            derivative_match_tol = 1e-12,
            derivative_match_strict = true)
    end

end # @testset Dispatcher

# DP.5.1  Mutation-proof procedure (verified 2026-05-13 before commit;
# see Phase-10 PathNetwork + Phase-11 BVP test files for the lineage):
#
#   Mutation A  --  in `dispatch_solve`'s BVPSegment branch, drop the
#     barycentric-recovered derivative at the right endpoint when setting
#     up the next IVP's IC: `new_cur_up = zero(CT)` instead of `up_at_end`.
#     Verified bite: 1 fail on DP.3.1 line 163 — third IVP segment starts
#     with u'=0 (mutated) instead of sinh(0.6)≈0.6367; u(0.75)≈1.199 vs
#     cosh(0.75)≈1.295, atol=1e-5 fails.  Both junctions' diagnostics
#     remain near-zero (the BVP terminus is correctly recovered; only the
#     downstream IVP IC is corrupted).  This isolates the BVP→IVP transition.
#
#   Mutation B  --  in `dispatch_solve`'s BVPSegment branch, swap u_a and
#     u_b in the bvp_solve call: pass `(seg.u_b, cur_u)` instead of
#     `(cur_u, seg.u_b)` as the 5th and 6th positional arguments.  The
#     BVP solves with swapped Dirichlet BCs; solution diverges from
#     cosh(t) over the whole interval.
#     Verified bite: 5 fails — DP.2.1 lines 116-117 (Δu = 0.140 since BC
#     no longer matches IVP terminus; Δu' = 0.941 vs ≤1e-7); DP.3.1 lines
#     153-154 (first junction Δu/Δu' identical); DP.3.1 line 163
#     (downstream cosh(0.75) mismatch — 1.011 vs 1.295).  This is the
#     analogue of Phase-11 BVP Mutation B (swap u_a/u_b in BC enforcement);
#     load-bearing across every test that exercises IVP→BVP coupling.
#
#   Mutation C  --  in `dispatch_solve`'s IVPSegment branch's dense-grid
#     append loop, gate on the first segment only: add `k == 1 || continue`
#     so segments k ≥ 2's IVP dense_grid values are silently dropped.
#     Verified bite: 1 fail + 2 errors on DP.3.1 — `idx_075 = findfirst(...)`
#     evaluates to `nothing` (third IVP's z=0.75 not in grid_z); the
#     subsequent `sol.grid_u[idx_075]` and `sol.grid_region[idx_075]`
#     indexing operations error with `MethodError` and `BoundsError`.
#     Confirms grid-concatenation across all IVP segments is load-bearing
#     for DP.3.1's downstream lookups.
#
# Restoration: all three mutations restored before commit per CLAUDE.md
# Rule 4.  Matches the Phase-6 + Phase-10 + Phase-11 inline pattern.
