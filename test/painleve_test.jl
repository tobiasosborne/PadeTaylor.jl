# test/painleve_test.jl -- ADR-0006 / bead `padetaylor-avl`.
#
# The PainleveProblem builder layer.  Tests fall in three groups:
#   PV.1.* -- constructor correctness: the assembled PadeTaylorProblem
#             carries the right RHS, frame, params, and (for the
#             transformed equations) ζ-frame IC.
#   PV.2.* -- fail-loud guards (CLAUDE.md Rule 1).
#   PV.3.* -- forwarding methods.
#   PV.4.* -- end-to-end smoke test through the path network.
#
# The RHS closures are exercised directly with plain-number arguments
# (they are also called with `Taylor1` args by the solver, but the
# closed-form value check is cleanest with scalars).
#
# Reference: docs/adr/0006-painleve-problem-layer.md; canonical forms
# FW2011 (PI), FW2014_*.md:47 (PII), ReegerFornberg2014_*.md:45 (PIV),
# CoordTransforms.jl / SheetTracker.jl docstrings (PIII/PV/PVI).

using Test
using PadeTaylor

@testset "Painleve (ADR-0006): PainleveProblem builder" begin

    @testset "PV.1.1: PI constructor -- RHS, frame, params, defaults" begin
        pp = PainleveProblem(:I; u0 = 0.1, up0 = 0.2, zspan = (0.0, 1.0))
        @test pp.equation == :I
        @test pp.frame == :direct
        @test pp.params == (;)
        @test pp.problem isa PadeTaylorProblem
        @test pp.problem.y0 == (0.1, 0.2)
        @test pp.problem.order == 30                       # default
        # u'' = 6u² + z.  At (z,u,up) = (2,3,1): 6·9 + 2 = 56.
        @test pp.problem.f(2.0, 3.0, 1.0) == 6 * 3.0^2 + 2.0
    end

    @testset "PV.1.2: PII constructor -- RHS carries α, order override" begin
        pp = PainleveProblem(:II; α = 0.7, u0 = 0.0, up0 = 1.0,
                             zspan = (0.0, 5.0), order = 20)
        @test pp.equation == :II
        @test pp.frame == :direct
        @test pp.params == (; α = 0.7)
        @test pp.problem.order == 20
        # u'' = 2u³ + zu + α.  At (z,u,up) = (2,3,1): 2·27 + 6 + 0.7.
        @test pp.problem.f(2.0, 3.0, 1.0) ≈ 2 * 3.0^3 + 2.0 * 3.0 + 0.7
        # α is genuinely wired in: f(0,0,0) == α.
        @test pp.problem.f(0.0, 0.0, 0.0) == 0.7
    end

    @testset "PV.1.3: PIV constructor -- RHS carries α,β; u=0 fails loud" begin
        pp = PainleveProblem(:IV; α = 1.0, β = 2.0, u0 = 0.5, up0 = 0.0,
                             zspan = (0.0, 3.0))
        @test pp.equation == :IV
        @test pp.frame == :direct
        @test pp.params == (; α = 1.0, β = 2.0)
        # u'' = (u')²/(2u) + (3/2)u³ + 4zu² + 2(z²−α)u + β/u.
        z, u, up = 1.0, 2.0, 3.0
        expected = up^2 / (2u) + (3//2) * u^3 + 4z * u^2 +
                   2 * (z^2 - 1.0) * u + 2.0 / u
        @test pp.problem.f(z, u, up) ≈ expected
        # u = 0 is a movable zero: the RHS must fail loud, not return Inf.
        @test_throws DomainError pp.problem.f(1.0, 0.0, 1.0)
    end

    @testset "PV.1.4: PIII constructor -- transformed frame, ζ-mapped IC" begin
        α, β, γ, δ = 1.0, 0.5, 0.25, 0.1
        z0, u0, up0 = 2.0 + 0im, 0.3 + 0im, 0.4 + 0im
        pp = PainleveProblem(:III; α = α, β = β, γ = γ, δ = δ,
                             u0 = u0, up0 = up0, zspan = (z0, 4.0 + 0im))
        @test pp.equation == :III
        @test pp.frame == :transformed
        @test pp.params == (; α = α, β = β, γ = γ, δ = δ)
        # The IC is mapped into the ζ-frame; the underlying problem is
        # built there.
        ζ0, w0, wp0 = pIII_z_to_ζ(z0, u0, up0)
        @test pp.problem.y0 == (w0, wp0)
        @test pp.problem.zspan[1] == ζ0
        # The RHS is the CoordTransforms factory verbatim -- sample match.
        rhs_ref = pIII_transformed_rhs(α, β, γ, δ)
        @test pp.problem.f(0.3 + 0.1im, 0.7 + 0im, 0.2 + 0im) ==
              rhs_ref(0.3 + 0.1im, 0.7 + 0im, 0.2 + 0im)
        # from_frame ∘ to_frame is the identity (round-trip).
        @test all(pp.from_frame(pp.to_frame(z0, u0, up0)...) .≈ (z0, u0, up0))
    end

    @testset "PV.1.5: PV + PVI constructors -- transformed frame, RHS match" begin
        α, β, γ, δ = 1.0, 0.5, 0.25, 0.1
        ppV = PainleveProblem(:V; α = α, β = β, γ = γ, δ = δ,
                              u0 = 0.6 + 0im, up0 = 0.2 + 0im,
                              zspan = (1.5 + 0im, 3.0 + 0im))
        @test ppV.equation == :V
        @test ppV.frame == :transformed
        rhsV = pV_transformed_rhs(α, β, γ, δ)
        @test ppV.problem.f(0.4 + 0im, 0.7 + 0im, 0.3 + 0im) ==
              rhsV(0.4 + 0im, 0.7 + 0im, 0.3 + 0im)

        # PVI reuses PV's z = exp(ζ) transform (FFW 2017 §2.2).
        ppVI = PainleveProblem(:VI; α = α, β = β, γ = γ, δ = δ,
                               u0 = 0.6 + 0im, up0 = 0.2 + 0im,
                               zspan = (2.0 + 0im, 3.0 + 0im))
        @test ppVI.equation == :VI
        @test ppVI.frame == :transformed
        rhsVI = pVI_transformed_rhs(α, β, γ, δ)
        @test ppVI.problem.f(0.4 + 0im, 0.7 + 0im, 0.3 + 0im) ==
              rhsVI(0.4 + 0im, 0.7 + 0im, 0.3 + 0im)
    end

    @testset "PV.2.1: fail-loud guards (CLAUDE.md Rule 1)" begin
        # Unknown equation.
        @test_throws ArgumentError PainleveProblem(:VII;
            u0 = 0.0, up0 = 1.0, zspan = (0.0, 1.0))
        # Missing required parameter.
        @test_throws ArgumentError PainleveProblem(:II;
            u0 = 0.0, up0 = 1.0, zspan = (0.0, 1.0))       # no α
        @test_throws ArgumentError PainleveProblem(:I;
            u0 = 0.0, zspan = (0.0, 1.0))                  # no up0
        # Unexpected parameter the equation does not take.
        @test_throws ArgumentError PainleveProblem(:I;
            α = 1.0, u0 = 0.0, up0 = 1.0, zspan = (0.0, 1.0))
        # Transformed-frame IC sitting on a fixed branch point.
        @test_throws ArgumentError PainleveProblem(:III;
            α = 1.0, β = 0.0, γ = 0.0, δ = 0.0,
            u0 = 1.0 + 0im, up0 = 0.0 + 0im,
            zspan = (0.0 + 0im, 2.0 + 0im))                # z = 0 branch pt
        @test_throws ArgumentError PainleveProblem(:VI;
            α = 1.0, β = 0.0, γ = 0.0, δ = 0.0,
            u0 = 1.0 + 0im, up0 = 0.0 + 0im,
            zspan = (1.0 + 0im, 3.0 + 0im))                # z = 1 branch pt
    end

    @testset "PV.3.1: forwarding returns PainleveSolution (ADR-0007)" begin
        # PII α=0, u(0)=0, u'(0)=1 -- smooth near the origin.
        pp = PainleveProblem(:II; α = 0.0, u0 = 0.0, up0 = 1.0,
                             zspan = (0.0, 0.5))

        # :direct forwarding wraps the raw solve output in a
        # PainleveSolution whose `.raw` is identical to the direct call.
        grid = ComplexF64[0.2 + 0.0im, 0.1 + 0.2im, -0.2 + 0.1im]
        sol_pp  = path_network_solve(pp, grid; h = 0.5)
        sol_raw = path_network_solve(pp.problem, grid; h = 0.5)
        @test sol_pp isa PainleveSolution
        @test sol_pp.equation == :II && sol_pp.frame == :direct
        @test sol_pp.raw.grid_u == sol_raw.grid_u
        @test sol_pp.raw.grid_up == sol_raw.grid_up

        sp_pp  = solve_pade(pp; h_max = 0.5)
        sp_raw = solve_pade(pp.problem; h_max = 0.5)
        @test sp_pp isa PainleveSolution
        @test sp_pp.raw.y == sp_raw.y

        # :transformed problems (ADR-0007 supersedes ADR-0006 ref. #2):
        # path_network_solve now works -- it maps the z-frame grid into
        # the ζ-frame, solves, and returns a PainleveSolution.  solve_pade
        # still throws for a complex ζ-domain (real-axis stepping is
        # undefined there) -- now pointing the caller at path_network_solve.
        ppIII = PainleveProblem(:III; α = 1.0, β = 0.5, γ = 0.25, δ = 0.1,
                                u0 = 0.3 + 0im, up0 = 0.4 + 0im,
                                zspan = (2.0 + 0im, 4.0 + 0im))
        solIII = path_network_solve(ppIII, ComplexF64[2.5 + 0im]; h = 0.5)
        @test solIII isa PainleveSolution
        @test solIII.frame == :transformed
        @test_throws ArgumentError solve_pade(ppIII; h_max = 0.5)
    end

    @testset "PV.4.1: end-to-end -- PII pole field is finite" begin
        # PII with α = 0, u(0) = 0, u'(0) = 1: a real solution, smooth
        # near the origin.  The path network over a small grid must
        # return finite values at every cell.  The forwarding method
        # returns a PainleveSolution; the grid is read via grid_values.
        pp = PainleveProblem(:II; α = 0.0, u0 = 0.0, up0 = 1.0,
                             zspan = (0.0, 2.0))
        grid = ComplexF64[x + im * y for x in -0.6:0.3:0.6
                                     for y in -0.6:0.3:0.6]
        sol = path_network_solve(pp, grid; h = 0.5)
        @test sol isa PainleveSolution
        gz, gu, gup = grid_values(sol)
        @test all(isfinite, real.(gu))
        @test all(isfinite, imag.(gu))
        @test all(isfinite, real.(gup))
        @test length(sol.raw.visited_z) ≥ 1
        @test gz == grid                                   # input order preserved
    end

end # @testset Painleve

# MUTATION-PROOF (verified 2026-05-14, bead `padetaylor-avl`):
#
#   Mutation A -- in `src/Painleve.jl`, perturb the PII RHS factory
#     `_pII_rhs(α) = (z,u,up) -> 2*u^3 + z*u + α`  →  `... + 2*α`.
#   Bite: PV.1.2's `f(0,0,0) == 0.7` goes RED (returns 1.4) and the
#   `f(2,3,1) ≈ …` assertion goes RED.  Restored.
#
#   Mutation B -- perturb the PIV RHS factory's `4 * z * u^2` term to
#     `5 * z * u^2`.
#   Bite: PV.1.3's `f(z,u,up) ≈ expected` goes RED.  Restored.
#
#   Mutation C -- in `src/Painleve.jl`'s `solve_pade(pp::PainleveProblem)`,
#     invert the transformed-frame guard condition: change
#     `pp.frame === :transformed && eltype(pp.problem.zspan) <: Complex`
#     to `pp.frame === :direct && …`.
#   Bite: PV.3.1's `@test_throws ArgumentError solve_pade(ppIII, …)` goes
#   RED -- the complex-ζ transformed problem no longer hits the guard, so
#   solve_pade runs and errors on `<(::Complex, ::Complex)` instead of
#   the expected ArgumentError.  Restored.
#
#   (Forwarding behaviour beyond the guard -- the PainleveSolution wrap,
#   the :transformed grid mapping, the z-frame access surface -- is
#   mutation-proven in `test/painleve_solution_test.jl`.)
