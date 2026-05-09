"""
Phase 5 RED tests for `PadeTaylor.PadeStepper.pade_step!`
(per DESIGN.md §4 Phase 5).

This is the big integration: Coefficients ↔ RobustPade ↔ StepControl
composed into one Padé-Taylor step.  After GREEN, the inner-loop
algorithm is proven correct; remaining phases (6-9) are scaffolding.

Oracle source (three-source triangulation):
    external/probes/padestepper-oracle/capture.wl  — Mathematica:
        closed-form WeierstrassP + NDSolve at WorkingPrecision=50
    external/probes/padestepper-oracle/capture.py  — mpmath.odefun
        at 40 dps (Taylor-method ODE solver, fully independent path)
    external/probes/padestepper-oracle/verify.jl   — three-way
        agreement enforcer at tolerances 1e-13 .. 1e-11
Pinned values: test/_oracle_padestepper.jl (auto-generated).

Test plan (deviations from DESIGN.md called out in each block):
    5.1.1  one ℘-step (z=0 → 0.5) on u''=6u² from FW ICs.
           u(0.5), u'(0.5) match the ℘ closed-form to 1e-13 rel.
    5.1.2  composition: continue from 5.1.1's state to z=0.9 (h=0.4).
           [DEVIATION: DESIGN says z=1.0, but z=1.0 is a pole of the
           closed-form — `1/u(1.0) = 0` per oracle.  Test at z=0.9
           instead, which lands close to the pole but safely before it.]
    5.1.3  one PI-step (z=0 → 0.5) on u''=6u²+z from FW ICs.
           u(0.5), u'(0.5) match the NDSolve oracle to 1e-13 rel.
           [DEVIATION: DESIGN says "tritronquée IC" without a paper-
           pinned value to high precision; we use FW ICs to test the
           same numerical pipeline + the +z RHS contribution.]
    5.1.4  near-pole step (z=0.9 → 0.95, h=0.05) on u''=6u² from
           u(0.9), u'(0.9) [pole at z=1].  [DEVIATION: DESIGN says
           "z=1.36" for the pole, but the actual pole on this
           lattice is at z = -c₁ = 1; `1/u(1) = 0` confirms.]
    5.1.5  Mutation-proof procedure documented at end (verified
           manually after GREEN).
"""

using Test
using PadeTaylor.PadeStepper: PadeStepperState, pade_step!

include("_oracle_padestepper.jl")

# Painlevé I right-hand side.
fPI(z, u, up) = 6 * u^2 + z

# Equianharmonic-Weierstrass right-hand side.
fW(z, u, up) = 6 * u^2

@testset "PadeStepper (Phase 5)" begin

    @testset "5.1.1 ℘ step h=0.5 from FW ICs (z=0 → 0.5)" begin
        # FW 2011-pinned ICs (FW2011_*.md:292-295) and the
        # three-source-pinned ℘ oracle at z=0.5.
        state = PadeStepperState{Float64}(0.0, u_0_FW, up_0_FW)
        pade_step!(state, fW, 30, 0.5)
        @test state.z ≈ 0.5
        @test isapprox(state.u,  u_5_1_1_at_05;  rtol = 1e-13)
        @test isapprox(state.up, up_5_1_1_at_05; rtol = 1e-13)
    end

    @testset "5.1.2 ℘ step composition (z=0 → 0.5 → 0.9)" begin
        # DEVIATION FROM DESIGN: target z=0.9 (h=0.4) rather than
        # z=1.0 (h=0.5).  z=1.0 is a pole on this lattice (verified
        # by oracle: inv_u_at_z1 == 0).  z=0.9 lands close to the
        # pole and exercises step composition + large-coefficient
        # behaviour.  rtol=1e-12 absorbs accumulated error from two
        # successive Padé-Taylor steps.
        state = PadeStepperState{Float64}(0.0, u_0_FW, up_0_FW)
        pade_step!(state, fW, 30, 0.5)   # → z=0.5
        pade_step!(state, fW, 30, 0.4)   # → z=0.9 (avoid pole at z=1)
        @test state.z ≈ 0.9
        @test isapprox(state.u,  u_5_1_2_at_09;  rtol = 1e-12)
        @test isapprox(state.up, up_5_1_2_at_09; rtol = 1e-12)
    end

    @testset "5.1.3 PI step h=0.5 from FW ICs (z=0 → 0.5)" begin
        # PI = u'' = 6u² + z.  Same ICs as 5.1.1; the +z RHS produces
        # a measurable contribution: u(0.5)_PI ≈ 4.034 vs u(0.5)_℘ ≈
        # 4.004; difference ≈ 0.0296 (oracle: diff_PI_minus_weier_at_05).
        state = PadeStepperState{Float64}(0.0, u_0_FW, up_0_FW)
        pade_step!(state, fPI, 30, 0.5)
        @test state.z ≈ 0.5
        @test isfinite(state.u) && isfinite(state.up)
        @test isapprox(state.u,  u_5_1_3_PI_at_05;  rtol = 1e-13)
        @test isapprox(state.up, up_5_1_3_PI_at_05; rtol = 1e-13)
        # Sanity: the +z RHS DID contribute (otherwise we'd just have ℘).
        @test abs(state.u  - u_5_1_1_at_05)  > 0.01
        @test abs(state.up - up_5_1_1_at_05) > 0.1
    end

    @testset "5.1.4 ℘ near-pole step (z=0.9 → 0.95, h=0.05)" begin
        # Pole at z=1 (since c₁ = -1 means u(z) = ℘(z + c₁) and ℘'s
        # principal pole is at the lattice origin).  Starting state at
        # z=0.9 has |u| ≈ 100, |u'| ≈ 2000 — the FW algorithm should
        # gracefully Padé-handle these large Taylor coefficients.
        # DEVIATION: DESIGN says "z ≈ 1.36" pole.  The actual pole on
        # this lattice is at z = -c₁ = 1 (oracle confirms: 1/u(1)=0).
        state = PadeStepperState{Float64}(0.9, u_state_at_09, up_state_at_09)
        pade_step!(state, fW, 30, 0.05)
        @test state.z ≈ 0.95
        @test isfinite(state.u) && isfinite(state.up)
        # rtol=1e-11: looser than 5.1.1 because the magnitudes are
        # ~400 and ~16000; absolute error budget at Float64 epsilon
        # times the magnitude is ~1e-13 × 16000 ≈ 1.6e-10.  We assert
        # 1e-11 which is achievable with a clean Padé.
        @test isapprox(state.u,  u_5_1_4_at_095;  rtol = 1e-11)
        @test isapprox(state.up, up_5_1_4_at_095; rtol = 1e-11)
    end
end

#=
Phase 5 mutation-proof procedure (test 5.1.5) — VERIFIED 2026-05-09
──────────────────────────────────────────────────────────────────
Both load-bearing evaluation paths (u, u') were independently mutation-
tested before commit; both mutations RED'd the test suite.

  Mutation A — sign flip on the Padé evaluation point for u
  (`src/PadeStepper.jl`, `pade_step!`, line ~224):
      new_u = _evaluate_pade(P_u, one_T)  →
      new_u = _evaluate_pade(P_u, -one_T)
    VERIFIED to fail the `state.u` assertion in all four testsets
    (5.1.1, 5.1.2, 5.1.3, 5.1.4).  After the h^k rescaling, t=1
    corresponds to z+h; t=-1 corresponds to z-h, giving u(z-h) ≠
    u(z+h).  Confirms the u-evaluation point is genuinely tested.

  Mutation C — flip the chain-rule h division on the u' path
  (`src/PadeStepper.jl`, `pade_step!`, line ~225):
      new_up = _evaluate_pade_deriv(P_u, one_T) / h_T  →
      new_up = _evaluate_pade_deriv(P_u, one_T) * h_T
    VERIFIED to fail the `state.up` assertion in all four testsets.
    The chain rule for u(z+h') = P_u(h'/h) gives u'(z+h) = P_u'(1)/h;
    multiplying by h instead of dividing yields the wrong scale by
    h² = 0.25 (5.1.1) or 0.0025 (5.1.4).  Confirms the chain-rule
    derivative-of-rescaled-Padé path is genuinely tested on both
    well-conditioned and near-pole cases.

Per worklog 001 + CLAUDE.md Rule 4 + Rule 5: each load-bearing
recurrence got an independent mutation that REDs at least one test
before this commit.  Restored before commit; full suite GREEN.
=#
