# test/ext_commonsolve_test.jl -- Phase 7 / bead `padetaylor-2vz` tests.
#
# Validates the `PadeTaylorCommonSolveExt` package extension per ADR-0003:
# `solve(prob, ::PadeTaylorAlg)` produces a `PadeTaylorSolution`
# bit-identical (modulo trivial evaluation-order differences) to a
# direct `solve_pade(prob; h_max, max_steps)` call on the same problem.
#
# Reference: docs/adr/0003-extensions-pattern.md (extensions pattern);
# ext/PadeTaylorCommonSolveExt.jl (the adapter being tested).
#
# v1 scope: 2nd-order PadeTaylorProblem only (matches solve_pade's v1
# branch).  1st-order branch errors with a Suggestion.
#
# Test groups:
#   CS.1.1  solve(prob, alg) ≡ solve_pade(prob)
#   CS.1.2  init + step! loop produces identical trajectory
#   CS.1.3  solve! produces same PadeTaylorSolution as solve
#   CS.2.1  Streaming: integrator's `done` flag flips at z_end
#   CS.2.2  Degenerate zspan: born-done integrator
#   CS.3.1  Fail-fast: h_max ≤ 0 throws ArgumentError
#   CS.3.2  Fail-fast: max_steps overrun throws ErrorException
#   CS.4.1  Mutation-proof commentary (verified manually)

using Test
using PadeTaylor
using CommonSolve

include(joinpath(@__DIR__, "_oracle_problems.jl"))

@testset "CommonSolveAdapter (Phase 7): SciML init/step!/solve!/solve" begin

    # Equianharmonic ℘ ODE: u'' = 6u^2, FW 2011 IC at z=0.  Single
    # pole-bridge segment with h_max=1.5 brackets the pole at z=1
    # (same setup as Phase-6 test 6.1.5's headline demo).
    fW(z, u, up) = 6 * u^2

    prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
    alg  = PadeTaylorAlg(; h_max = 1.5)

    @testset "CS.1.1: solve(prob, alg) ≡ solve_pade(prob; h_max=alg.h_max)" begin
        sol_cs    = CommonSolve.solve(prob, alg)
        sol_pade  = solve_pade(prob; h_max = 1.5)
        # Parallel-vector equality.
        @test sol_cs.z    == sol_pade.z
        @test sol_cs.y    == sol_pade.y
        @test sol_cs.h    == sol_pade.h
        # Padé approximants: compare via the dense-eval callable at sample t.
        @test length(sol_cs.pade) == length(sol_pade.pade)
        for t_eval in (0.3, 0.5, 0.7)
            @test sol_cs(t_eval) == sol_pade(t_eval)
        end
    end

    @testset "CS.1.2: init + step! loop ≡ direct solve" begin
        # Build the same trajectory by hand via the streaming API.
        integ = CommonSolve.init(prob, alg)
        @test integ.steps == 0
        @test !integ.done
        @test length(integ.z_vec) == 1   # IC pre-pushed
        @test length(integ.h_vec) == 0
        n_loop = 0
        while !integ.done
            n_loop += 1
            CommonSolve.step!(integ)
            n_loop > 1000 && error("test loop runaway")
        end
        @test integ.done
        # The trajectory matches solve_pade.
        sol_pade = solve_pade(prob; h_max = 1.5)
        @test integ.z_vec == sol_pade.z
        @test integ.y_vec == sol_pade.y
        @test integ.h_vec == sol_pade.h
        @test length(integ.pade_vec) == length(sol_pade.pade)
    end

    @testset "CS.1.3: solve! after manual init ≡ solve" begin
        integ    = CommonSolve.init(prob, alg)
        sol_via_solve!  = CommonSolve.solve!(integ)
        sol_via_solve   = CommonSolve.solve(prob, alg)
        sol_pade        = solve_pade(prob; h_max = 1.5)
        @test sol_via_solve!.z == sol_pade.z
        @test sol_via_solve!.y == sol_pade.y
        @test sol_via_solve.z  == sol_pade.z
        @test sol_via_solve.y  == sol_pade.y
        # Dense-eval sample point.
        @test sol_via_solve!(1.0) == sol_pade(1.0)
        @test sol_via_solve(1.0)  == sol_pade(1.0)
    end

    @testset "CS.2.1: streaming — done flag flips at z_end, uneven final step" begin
        # Uneven final step exercise: zspan (0, 0.55) with h_max=0.2.
        # Even loop: 0.2 + 0.2 + 0.15 = 0.55.  The MIN clamp on the final
        # step is load-bearing — without it the integrator would overshoot
        # to z=0.6.  This is the mutation-G target.
        prob_uneven = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 0.55);
                                        order = 30)
        alg_uneven  = PadeTaylorAlg(; h_max = 0.2)
        integ       = CommonSolve.init(prob_uneven, alg_uneven)
        @test !integ.done
        steps_done = 0
        while !integ.done
            steps_done += 1
            CommonSolve.step!(integ)
            steps_done > 100 && error("test runaway")
        end
        @test integ.done
        @test steps_done == 3
        @test integ.z_vec[end] ≈ 0.55           # exactly z_end, not 0.6
        @test integ.h_vec[end] ≈ 0.15           # final step clamped from 0.2
        # step! after done is a no-op.
        prev_z_len = length(integ.z_vec)
        CommonSolve.step!(integ)
        @test length(integ.z_vec) == prev_z_len
    end

    @testset "CS.2.2: degenerate zspan → born-done integrator" begin
        # zspan with z_start == z_end -> the constructor rejects this in
        # PadeTaylorProblem, so we use the smallest practical span.  But
        # the integrator should still handle the case where z_start ≥ z_end
        # (e.g., a user-supplied prob with reversed bounds).  Since the
        # constructor enforces z_start != z_end, we test the constructor's
        # rejection here instead — the integrator's born-done branch is
        # unreachable through public API but kept for defence-in-depth.
        @test_throws ArgumentError PadeTaylorProblem(fW, (u_0_FW, up_0_FW),
                                                    (0.5, 0.5); order = 30)
    end

    @testset "CS.3.1: fail-fast — h_max ≤ 0 throws at init" begin
        # PadeTaylorAlg constructor doesn't validate; init does.
        bad_alg_neg = PadeTaylorAlg(; h_max = -0.5)
        @test_throws ArgumentError CommonSolve.init(prob, bad_alg_neg)
        bad_alg_zero = PadeTaylorAlg(; h_max = 0.0)
        @test_throws ArgumentError CommonSolve.init(prob, bad_alg_zero)
    end

    @testset "CS.3.2: fail-fast — max_steps overrun throws at step!" begin
        # max_steps too small for the segment ⇒ step! errors mid-loop.
        # zspan (0, 1.5) with h_max=0.1 → ~15 steps; max_steps=3 ⇒ overrun.
        prob_long  = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5);
                                       order = 30)
        alg_short  = PadeTaylorAlg(; h_max = 0.1, max_steps = 3)
        integ      = CommonSolve.init(prob_long, alg_short)
        # First 3 step!s succeed; 4th throws.
        CommonSolve.step!(integ)
        CommonSolve.step!(integ)
        CommonSolve.step!(integ)
        @test integ.steps == 3
        @test !integ.done
        @test_throws ErrorException CommonSolve.step!(integ)
    end

end # @testset CommonSolveAdapter

# CS.4.1  Mutation-proof procedure (verified 2026-05-13):
#
#   Mutation G  --  in `step!`, change `h_step = min(h_max_T, z_end - state.z)`
#     to `h_step = h_max_T` (never clamp the final step to land on z_end).
#     Verified bite: 2 fails on CS.2.1 — `integ.z_vec[end] ≈ 0.55` evaluates
#     `0.6 ≈ 0.55` (overshoot); `integ.h_vec[end] ≈ 0.15` evaluates `0.2 ≈ 0.15`
#     (unclamped step length).  The min-clamp on the final step is what
#     ensures the integrator lands EXACTLY on `zspan[2]`.
#
#   Mutation H  --  in `step!`, drop the `integ.done && return integ`
#     early-out so step! attempts another step after done.  The inner
#     `pade_step_with_pade!` accepts h_step=0 (since z_end - state.z = 0
#     at that point) and silently appends a zero-length segment.
#     Verified bite: 1 fail on CS.2.1 — `length(integ.z_vec)` after the
#     extra step! is 5 instead of 4 (the IC + 3 real steps); the no-op-on-
#     done contract is broken.
#
# Restoration: both mutations restored before commit per CLAUDE.md Rule 4.
