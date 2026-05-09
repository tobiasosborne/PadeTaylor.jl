using Test
using PadeTaylor

@testset "PadeTaylor.jl" begin
    @testset "umbrella loads" begin
        # Confirms the umbrella module + its 6 sub-modules load cleanly.
        @test isdefined(PadeTaylor, :LinAlg)
        @test isdefined(PadeTaylor, :RobustPade)
        @test isdefined(PadeTaylor, :Coefficients)
        @test isdefined(PadeTaylor, :StepControl)
        @test isdefined(PadeTaylor, :PadeStepper)
        @test isdefined(PadeTaylor, :Problems)
    end

    # Per-module test files. Each file is a self-contained @testset.
    # Phases 1–6 fill in these test files in turn.
    include("linalg_test.jl")
    include("robustpade_test.jl")
    include("coefficients_test.jl")
    include("stepcontrol_test.jl")
    include("padestepper_test.jl")
    include("problems_test.jl")
end
