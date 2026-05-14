# test/painleve_solution_test.jl -- ADR-0007 / bead `padetaylor-p1b`.
#
# The PainleveSolution wrapper: a self-describing solve-output object
# carrying the Painlevé identity (equation, params, named-solution tag,
# coordinate frame) around whichever raw type the solver produced, with
# a uniform z-frame access surface (callable, poles, grid_values).
#
# Test groups:
#   PS.1.* -- provenance fields + accessors.
#   PS.2.* -- the callable: :direct transparent forwarding, :transformed
#             frame round-trip, fail-loud on a grid-type raw.
#   PS.3.* -- poles(): :direct forwarding identity, :transformed mapping.
#   PS.4.* -- show: multi-line provenance + one-line form.
#   PS.5.* -- grid_values(): :direct passthrough, :transformed mapping.
#   PS.6.* -- fail-loud on an unwired raw type (Rule 1).
#
# Several PS.3/PS.5 tests wrap a *hand-built* PathNetworkSolution with a
# known single pole and known grid point: the wrapper's job is to
# forward / frame-map, so a deterministic raw makes the mapping
# invariants exact and non-vacuous (the same idiom as polefield PF.2.1).
#
# Reference: docs/adr/0007-painleve-solution-wrapper.md.

using Test
using PadeTaylor

@testset "PainleveSolution (ADR-0007): self-describing solve-output wrapper" begin

    # A hand-built single-node PathNetworkSolution whose stored Padé has
    # denominator D(t) = 1 - 0.5t -> a simple pole at t* = 2, numerator
    # N(t) = 1 -> residue N(2)/D'(2) = -2 (not a Froissart doublet).  The
    # node sits at z_v = 0.5 with canonical step h = 0.5, so extract_poles
    # surfaces one pole at z = 0.5 + 0.5·2 = 1.5.  The grid carries one
    # point (z, u, u') = (0.3, 1.2, 0.7).
    T  = Float64
    Pp = PadeApproximant{Complex{T}}(Complex{T}[1.0], Complex{T}[1.0, -0.5],
                                     0, 1, one(Complex{T}))
    pns = PathNetworkSolution{T}(
        Complex{T}[0.5],          # visited_z
        Complex{T}[1.0],          # visited_u
        Complex{T}[0.0],          # visited_up
        [Pp],                     # visited_pade
        T[0.5],                   # visited_h
        Int[0],                   # visited_parent
        Complex{T}[0.3],          # grid_z
        Complex{T}[1.2],          # grid_u
        Complex{T}[0.7],          # grid_up
    )

    @testset "PS.1.1: provenance fields + accessors" begin
        pp  = PainleveProblem(:I; u0 = 0.1, up0 = 0.2, zspan = (0.0, 1.0))
        sol = solve_pade(pp; h_max = 0.5)
        @test sol isa PainleveSolution
        @test sol.raw isa PadeTaylorSolution
        @test equation(sol)     == :I
        @test parameters(sol)   == (;)
        @test solutionname(sol) === nothing      # generic IC, not a named branch
        @test sol.frame         == :direct

        ppII = PainleveProblem(:II; α = 0.7, u0 = 0.0, up0 = 1.0,
                               zspan = (0.0, 0.5))
        solII = solve_pade(ppII; h_max = 0.5)
        @test equation(solII)   == :II
        @test parameters(solII) == (; α = 0.7)
    end

    @testset "PS.2.1: :direct callable forwards transparently" begin
        # PII α=0, u(0)=0, u'(0)=1 -- smooth near the origin.
        pp  = PainleveProblem(:II; α = 0.0, u0 = 0.0, up0 = 1.0,
                              zspan = (0.0, 0.6))
        sol = solve_pade(pp; h_max = 0.3)
        # The wrapper callable is byte-identical to the raw callable for
        # a :direct problem.
        for z in (0.0, 0.15, 0.42, 0.6)
            @test sol(z) == sol.raw(z)
        end
        # And it reproduces the initial condition at z0 (the IC point).
        u0, up0 = sol(0.0)
        @test u0  ≈ 0.0 atol = 1e-12
        @test up0 ≈ 1.0 atol = 1e-12
    end

    @testset "PS.2.2: grid-type raw is not callable (fail loud, Rule 1)" begin
        # A PathNetworkSolution is a grid scatter, not a dense
        # interpolant -- the wrapper callable must throw, not fake a value.
        pp  = PainleveProblem(:II; α = 0.0, u0 = 0.0, up0 = 1.0,
                              zspan = (0.0, 0.5))
        sol = path_network_solve(pp, ComplexF64[0.2 + 0.0im]; h = 0.5)
        @test sol.raw isa PathNetworkSolution
        @test_throws ArgumentError sol(0.2)
    end

    @testset "PS.2.3: :transformed callable round-trips through the ζ-frame" begin
        # PV with real ICs/span -> a real ζ-domain, so solve_pade applies.
        # The callable maps z -> ζ, evaluates the ζ-frame trajectory, and
        # maps (w,w') -> (u,u') back.  At the IC point that whole pipeline
        # must reproduce the initial condition -- a genuine, non-circular
        # round-trip check across the constructor's z->ζ map, the solve,
        # and the callable's ζ->z map.
        ppV = PainleveProblem(:V; α = 0.1, β = 0.1, γ = 0.1, δ = 0.1,
                              u0 = 0.5, up0 = 0.1, zspan = (2.0, 3.0))
        sol = solve_pade(ppV; h_max = 0.2)
        @test sol isa PainleveSolution
        @test sol.frame == :transformed
        u0, up0 = sol(2.0)                      # the IC point z0 = 2.0
        @test u0  ≈ 0.5 atol = 1e-9
        @test up0 ≈ 0.1 atol = 1e-9
        # An interior z evaluates finite (and stays real for a real solve).
        ui, upi = sol(2.5)
        @test isfinite(ui) && isfinite(upi)
    end

    @testset "PS.3.1: poles() forwards identically for a :direct problem" begin
        # :direct frame -> from_frame is the identity, so poles(sol) must
        # be exactly extract_poles(sol.raw).  Pinned on the hand-built
        # single-pole network so the result is known and non-empty.
        sol = PainleveSolution(:I, (;), nothing, :direct, tuple, tuple, pns)
        p_wrap = poles(sol; min_support = 1)
        p_raw  = extract_poles(pns; min_support = 1)
        @test !isempty(p_wrap)
        @test p_wrap == p_raw
        @test p_wrap[1] ≈ 1.5 + 0.0im            # z_v + h·t* = 0.5 + 0.5·2
    end

    @testset "PS.3.2: poles() maps back to z-frame for a :transformed problem" begin
        # Same raw network, wrapped as a :transformed PV problem.  PV's
        # inverse coordinate map is z = exp(ζ), so each ζ-frame pole must
        # come back as exp(ζ).
        dir   = PainleveSolution(:V, (;), nothing, :direct,
                                 tuple, tuple, pns)
        trans = PainleveSolution(:V, (;), nothing, :transformed,
                                 pV_z_to_ζ, pV_ζ_to_z, pns)
        p_dir   = poles(dir;   min_support = 1)
        p_trans = poles(trans; min_support = 1)
        @test !isempty(p_dir)
        @test p_trans ≈ exp.(p_dir)              # z = exp(ζ), the PV map
    end

    @testset "PS.4.1: show -- multi-line provenance + one-line form" begin
        pp  = PainleveProblem(:II; α = 0.7, u0 = 0.0, up0 = 1.0,
                              zspan = (0.0, 0.5))
        sol = solve_pade(pp; h_max = 0.5)

        full = sprint(show, MIME("text/plain"), sol)
        @test occursin("Painlevé II", full)
        @test occursin("α = 0.7", full)
        @test occursin("direct", full)
        @test occursin("generic IC", full)       # name === nothing
        @test occursin("trajectory", full)        # raw summary

        line = sprint(show, sol)
        @test occursin("PainleveSolution", line)
        @test occursin("Painlevé II", line)
    end

    @testset "PS.5.1: grid_values -- :direct passthrough" begin
        # :direct frame -> grid_values returns the raw point arrays as-is.
        sol = PainleveSolution(:I, (;), nothing, :direct, tuple, tuple, pns)
        gz, gu, gup = grid_values(sol)
        @test gz  == pns.grid_z
        @test gu  == pns.grid_u
        @test gup == pns.grid_up
    end

    @testset "PS.5.2: grid_values -- :transformed maps point values back" begin
        # Wrapped as :transformed PV.  pV_ζ_to_z(ζ,w,w') = (exp ζ, w, w'/z),
        # so the z-frame grid is exp(ζ-grid), u is unchanged (u = w), and
        # u' is rescaled by 1/z.
        sol = PainleveSolution(:V, (;), nothing, :transformed,
                               pV_z_to_ζ, pV_ζ_to_z, pns)
        gz, gu, gup = grid_values(sol)
        @test gz  ≈ exp.(pns.grid_z)
        @test gu  ≈ pns.grid_u
        @test gup ≈ pns.grid_up ./ exp.(pns.grid_z)
    end

    @testset "PS.5.3: :transformed path_network_solve round-trips the grid" begin
        # End-to-end through the forwarding method (not a hand-built raw):
        # the caller's z-frame grid is mapped into the ζ-frame to solve,
        # and grid_values must hand the *same z-frame grid* back.  This is
        # the load-bearing invariant for the :transformed grid mapping in
        # `path_network_solve(::PainleveProblem, …)`.
        ppV = PainleveProblem(:V; α = 0.1, β = 0.1, γ = 0.1, δ = 0.1,
                              u0 = 0.5 + 0.0im, up0 = 0.1 + 0.0im,
                              zspan = (2.0 + 0.0im, 3.0 + 0.0im))
        grid = ComplexF64[2.2 + 0.1im, 2.5 - 0.1im, 2.7 + 0.0im]
        sol  = path_network_solve(ppV, grid; h = 0.4)
        @test sol isa PainleveSolution && sol.frame == :transformed
        gz, gu, gup = grid_values(sol)
        @test gz ≈ grid                          # z-frame grid round-trips
        @test all(isfinite, gu) && all(isfinite, gup)
    end

    @testset "PS.6.1: unwired raw type fails loud on every accessor (Rule 1)" begin
        # A PainleveSolution whose raw is neither callable nor a known
        # pole/grid source: each access verb must throw with a suggestion,
        # never silently return a wrong-but-plausible value.
        bogus = PainleveSolution(:I, (;), nothing, :direct,
                                 tuple, tuple, "not a solution")
        @test_throws ArgumentError bogus(0.5)         # callable
        @test_throws ArgumentError poles(bogus)       # poles
        @test_throws ArgumentError grid_values(bogus) # grid_values
    end

end # @testset PainleveSolution

# ----------------------------------------------------------------------
# Mutation-proof procedure (CLAUDE.md Rule 3 / Rule 4).
#
# Each mutation below was applied to `src/PainleveSolution.jl` or
# `src/Painleve.jl`, the suite confirmed RED, then the mutation reverted.
#
#   M1 -- frame map projection.  In `src/PainleveSolution.jl`, change
#         `_coord(map, c) = map(c, zero(c), zero(c))[1]` to `[2]`.
#         Bite: PS.2.3 (the :transformed IC round-trip no longer
#         reproduces u0/up0), PS.3.2 (poles map through the wrong
#         component).  Restored.
#
#   M2 -- grid-type callable guard.  Delete the `_assert_callable(sol.raw)`
#         line from the `(sol::PainleveSolution)(z)` callable.
#         Bite: PS.2.2's `@test_throws ArgumentError sol(0.2)` goes RED
#         and PS.6.1's `@test_throws ArgumentError bogus(0.5)` goes RED --
#         both calls now fall through to `sol.raw(z)`, throwing a
#         MethodError (PathNetworkSolution / String not callable) rather
#         than the expected ArgumentError.  Restored.
#
#   M3 -- :transformed grid mapping in the forwarding method.  In
#         `src/Painleve.jl`'s `path_network_solve(pp::PainleveProblem, …)`,
#         in the `:transformed` branch use `grid` directly instead of the
#         mapped `ζgrid`.
#         Bite: PS.5.3's `gz ≈ grid` goes RED -- the solve runs on
#         z-frame points treated as ζ, so the round-trip
#         `grid_values` z-array becomes `exp.(grid)` instead of the
#         input grid.  (PS.5.2 uses a hand-built raw and so does *not*
#         exercise the forwarding method; PS.5.3 is the end-to-end test
#         that pins this mutation.)  Restored.
#
#   M4 -- poles() not-wired guard.  In `src/PainleveSolution.jl`'s
#         `poles`, change the `raw isa … || raw isa …` guard to `true`.
#         Bite: PS.6.1's `@test_throws ArgumentError poles(bogus)` goes
#         RED (it now calls `extract_poles("not a solution")` -> a
#         MethodError, not the expected ArgumentError).  Restored.
# ----------------------------------------------------------------------
