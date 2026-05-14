# test/painleve_named_test.jl -- ADR-0008 / bead `padetaylor-c86`.
#
# Named-transcendent constructors: tritronquee(:I) and hastings_mcleod().
# Each returns a PainleveProblem carrying a literature-pinned initial
# condition + a `name` tag that propagates into PainleveSolution.name.
#
# Test groups:
#   NT.1.* -- constructor correctness: the built PainleveProblem carries
#             the cited IC (asserted verbatim against the reference
#             values -- the IC *is* the oracle, CLAUDE.md Rule 5), the
#             right equation/params/frame, and the `name` tag.
#   NT.2.* -- end-to-end: `name` propagates through the solvers into
#             solutionname(sol); the Hastings-McLeod branches are
#             negatives of each other (the PII alpha=0 symmetry
#             u -> -u -- a genuine cross-check); HM decays on x > 0.
#   NT.3.* -- fail-loud guards (CLAUDE.md Rule 1).
#
# Ground truth (verified verbatim at coding time):
#   tritronquee  -- references/markdown/FW2011_painleve_methodology_JCP230/
#                   FW2011_painleve_methodology_JCP230.md:224-229
#   hastings_mcleod -- references/markdown/
#                   FW2014_second_PII_exploration_FoCM14/
#                   FW2014_second_PII_exploration_FoCM14.md:252-258
#
# Reference: docs/adr/0008-named-transcendent-constructors.md.

using Test
using PadeTaylor

@testset "PainleveNamed (ADR-0008): named transcendent constructors" begin

    @testset "NT.1.1: tritronquee(:I) -- FW 2011 eq. (4.1) IC + name tag" begin
        pp = tritronquee(:I)
        @test pp isa PainleveProblem
        @test pp.equation == :I
        @test pp.frame    == :direct
        @test pp.params   == (;)
        @test pp.name     == :tritronquee
        # The IC is the oracle: FW 2011 §4.1, verbatim to 16 digits.
        @test pp.problem.y0 == (-0.1875543083404949, 0.3049055602612289)
        @test pp.problem.order == 30                       # default
        # zspan + order are caller-chosen.
        pp2 = tritronquee(:I; zspan = (0.0, 25.0), order = 40)
        @test pp2.problem.order == 40
        @test pp2.problem.zspan == (0.0, 25.0)
        @test pp2.name == :tritronquee
    end

    @testset "NT.1.2: hastings_mcleod() -- FW 2014 IC, both branches" begin
        # :positive (default): the (+, -) sign pattern.
        pp = hastings_mcleod()
        @test pp isa PainleveProblem
        @test pp.equation == :II
        @test pp.frame    == :direct
        @test pp.params   == (; α = 0.0)
        @test pp.name     == :hastings_mcleod
        @test pp.problem.y0 == (0.3670615515480784, -0.2953721054475501)

        # :negative: the (-, +) sign pattern -- the sign-symmetric copy.
        ppn = hastings_mcleod(; branch = :negative)
        @test ppn.problem.y0 == (-0.3670615515480784, 0.2953721054475501)
        @test ppn.name == :hastings_mcleod
    end

    @testset "NT.2.1: `name` propagates through the solvers" begin
        # solve_pade -> PainleveSolution carries the tag.
        sol = solve_pade(tritronquee(:I; zspan = (0.0, 0.5)); h_max = 0.5)
        @test sol isa PainleveSolution
        @test solutionname(sol) == :tritronquee
        @test equation(sol) == :I

        # path_network_solve -> likewise.
        soln = path_network_solve(hastings_mcleod(),
                                  ComplexF64[0.2 + 0.0im]; h = 0.5)
        @test soln isa PainleveSolution
        @test solutionname(soln) == :hastings_mcleod
        @test equation(soln) == :II
    end

    @testset "NT.2.2: Hastings-McLeod branches are negatives (PII α=0 u→-u)" begin
        # PII at α = 0 is u'' = 2u³ + zu, invariant under u -> -u.  The
        # :positive and :negative HM solutions are therefore exact
        # negatives -- a genuine cross-check that the branch sign is
        # wired through both the IC *and* the integration.
        sp = solve_pade(hastings_mcleod(; branch = :positive,
                                        zspan = (0.0, 2.0)); h_max = 0.5)
        sn = solve_pade(hastings_mcleod(; branch = :negative,
                                        zspan = (0.0, 2.0)); h_max = 0.5)
        for z in (0.0, 0.4, 1.1, 1.8)
            up_pos, upp_pos = sp(z)
            up_neg, upp_neg = sn(z)
            @test up_pos  ≈ -up_neg  atol = 1e-10
            @test upp_pos ≈ -upp_neg atol = 1e-10
        end
    end

    @testset "NT.2.3: Hastings-McLeod decays on the positive real axis" begin
        # HM (α=0) is pole-free and non-oscillatory on ℝ, with
        # u(x) ~ Ai(x) -> 0 as x -> +∞ (FW 2014).  Integrating the
        # :positive branch from 0, |u| must shrink from its IC value.
        sol = solve_pade(hastings_mcleod(; zspan = (0.0, 3.0)); h_max = 0.5)
        u0, _ = sol(0.0)
        u3, _ = sol(3.0)
        @test u0 ≈ 0.3670615515480784 atol = 1e-10   # IC reproduced
        @test isfinite(u3)
        @test 0.0 < u3 < u0                          # decaying toward 0
    end

    @testset "NT.3.1: fail-loud guards (CLAUDE.md Rule 1)" begin
        # tritronquee only exists in-tree for :I.
        @test_throws ArgumentError tritronquee(:II)
        @test_throws ArgumentError tritronquee(:V)
        # IC point is fixed at z = 0; a mismatched zspan[1] must throw.
        @test_throws ArgumentError tritronquee(:I; zspan = (1.0, 10.0))
        @test_throws ArgumentError hastings_mcleod(; zspan = (2.0, 5.0))
        # Unknown Hastings-McLeod branch.
        @test_throws ArgumentError hastings_mcleod(; branch = :sideways)
    end

end # @testset PainleveNamed

# ----------------------------------------------------------------------
# Mutation-proof procedure (CLAUDE.md Rule 3 / Rule 4).
#
# Each mutation applied to `src/PainleveNamed.jl` or `src/PainleveSolution
# .jl`, the suite confirmed RED, then reverted.
#
#   M1 -- the tritronquée IC.  In `src/PainleveNamed.jl`, perturb the
#         last digit of `u0`: `-0.1875543083404949` -> `-0.1875543083404948`.
#         Bite: NT.1.1's verbatim `pp.problem.y0 == (...)` goes RED.
#         Restored.
#
#   M2 -- the `name` plumbing.  In `src/PainleveSolution.jl`, revert
#         `_painleve_solution`'s `name = pp.name` default back to
#         `name = nothing`.
#         Bite: NT.2.1's `solutionname(sol) == :tritronquee` and
#         `== :hastings_mcleod` both go RED (the tag never reaches the
#         PainleveSolution).  Restored.
#
#   M3 -- the Hastings-McLeod branch sign.  In `src/PainleveNamed.jl`,
#         make `:negative` not flip the sign: `s = branch === :positive
#         ? 1.0 : -1.0` -> `s = 1.0`.
#         Bite: NT.1.2's `:negative` `ppn.problem.y0 == (...)` goes RED,
#         and NT.2.2's branch-symmetry assertions go RED (the two
#         "branches" are now identical, not negatives).  Restored.
# ----------------------------------------------------------------------
