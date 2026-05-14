# test/ext_makie_test.jl -- bead `padetaylor-ylr` tests.
#
# Validates the `PadeTaylorMakieExt` package extension per ADR-0003 +
# ADR-0007: `painleveplot(::PainleveSolution)` builds a complex-z-plane
# Makie figure (trajectory / grid + poles) from the wrapper's public
# surface only.  The extension is presentation-only, so the tests
# inspect the *structure* of the returned figure -- a Figure with one
# titled Axis, the expected number of plot objects -- rather than
# rendering pixels (no backend is loaded; Makie builds the scene graph
# regardless).
#
# Test groups:
#   MK.1.1  trajectory-backed PainleveSolution -> titled Figure + Axis.
#   MK.1.2  grid-backed PainleveSolution -> Figure.
#   MK.1.3  show_poles toggles the pole overlay (pinned on a hand-built
#           raw with a known pole so the overlay is non-vacuous).
#
# Reference: ext/PadeTaylorMakieExt.jl;
# docs/adr/0007-painleve-solution-wrapper.md.

using Test
using PadeTaylor
using Makie

@testset "MakieExt (bead padetaylor-ylr): painleveplot recipe" begin

    # Hand-built single-node PathNetworkSolution with one known pole at
    # z = 0.5 + 0.5·2 = 1.5 (D(t) = 1 - 0.5t -> root t* = 2; the same
    # fixture polefield PF.2.1 / painleve_solution PS.3.* use).
    T  = Float64
    Pp = PadeApproximant{Complex{T}}(Complex{T}[1.0], Complex{T}[1.0, -0.5],
                                     0, 1, one(Complex{T}))
    pns = PathNetworkSolution{T}(
        Complex{T}[0.5], Complex{T}[1.0], Complex{T}[0.0],
        [Pp], T[0.5], Int[0],
        Complex{T}[0.3, 0.6], Complex{T}[1.2, 1.4], Complex{T}[0.7, 0.9],
    )

    @testset "MK.1.1: trajectory-backed solution → titled Figure + Axis" begin
        pp  = PainleveProblem(:II; α = 0.0, u0 = 0.0, up0 = 1.0,
                              zspan = (0.0, 0.6))
        sol = solve_pade(pp; h_max = 0.3)
        fig = painleveplot(sol)
        @test fig isa Makie.Figure
        ax = fig.content[1]
        @test ax isa Makie.Axis
        @test occursin("Painlevé II", ax.title[])
        # PadeTaylorSolution raw draws a trajectory line + an IC marker.
        @test length(ax.scene.plots) ≥ 2
    end

    @testset "MK.1.2: grid-backed solution → Figure" begin
        pp   = PainleveProblem(:II; α = 0.0, u0 = 0.0, up0 = 1.0,
                               zspan = (0.0, 0.5))
        grid = ComplexF64[x + im * y for x in -0.4:0.2:0.4
                                     for y in -0.4:0.2:0.4]
        sol  = path_network_solve(pp, grid; h = 0.5)
        fig  = painleveplot(sol)
        @test fig isa Makie.Figure
        @test fig.content[1] isa Makie.Axis
    end

    @testset "MK.1.3: show_poles toggles the pole overlay" begin
        # Wrap the hand-built single-pole network.  With min_support = 1
        # the pole at z = 1.5 is found, so show_poles = true adds one
        # extra scatter over show_poles = false.
        sol = PainleveSolution(:I, (;), nothing, :direct, tuple, tuple, pns)

        fig_on  = painleveplot(sol; show_poles = true,
                               pole_kwargs = (; min_support = 1))
        fig_off = painleveplot(sol; show_poles = false)

        n_on  = length(fig_on.content[1].scene.plots)
        n_off = length(fig_off.content[1].scene.plots)
        @test n_on == n_off + 1            # the pole scatter is the only delta
    end

end # @testset MakieExt

# ----------------------------------------------------------------------
# Mutation-proof procedure (CLAUDE.md Rule 3 / Rule 4).
#
# Applied to `ext/PadeTaylorMakieExt.jl`, suite confirmed RED, reverted.
#
#   M1 -- pole overlay.  In `painleveplot`, change `if show_poles` to
#         `if false`.  Bite: MK.1.3's `n_on == n_off + 1` goes RED
#         (`n_on == n_off` -- the pole scatter is never drawn).  Restored.
#
#   M2 -- raw-type dispatch.  Delete the `_draw_raw!(ax, ::PadeTaylor
#         Solution, …)` method so the trajectory falls through to no
#         drawing.  Bite: MK.1.1's `length(ax.scene.plots) ≥ 2` goes RED
#         (a MethodError, actually -- the call has no applicable method),
#         confirming the trajectory branch is load-bearing.  Restored.
# ----------------------------------------------------------------------
