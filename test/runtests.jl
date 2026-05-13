using Test
using PadeTaylor

@testset "PadeTaylor.jl" begin
    @testset "umbrella loads" begin
        # Confirms the umbrella module + its 8 sub-modules load cleanly.
        @test isdefined(PadeTaylor, :LinAlg)
        @test isdefined(PadeTaylor, :RobustPade)
        @test isdefined(PadeTaylor, :Coefficients)
        @test isdefined(PadeTaylor, :StepControl)
        @test isdefined(PadeTaylor, :PadeStepper)
        @test isdefined(PadeTaylor, :Problems)
        @test isdefined(PadeTaylor, :PathNetwork)
        @test isdefined(PadeTaylor, :BVP)
        @test isdefined(PadeTaylor, :Dispatcher)
    end

    # Per-module test files. Each file is a self-contained @testset.
    # Phases 1–6 + Phase 10 PathNetwork + Phase 11 BVP + Phase 12
    # Dispatcher fill in these.
    include("linalg_test.jl")
    include("robustpade_test.jl")
    include("coefficients_test.jl")
    include("stepcontrol_test.jl")
    include("padestepper_test.jl")
    include("problems_test.jl")
    include("pathnetwork_test.jl")
    include("bvp_test.jl")
    include("dispatcher_test.jl")
    include("ext_commonsolve_test.jl")
    include("ext_arblib_test.jl")
end
