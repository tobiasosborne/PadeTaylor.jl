"""
V3c tests for `PadeTaylor.VectorStepControl` (bead `padetaylor-0ln.10`;
v0.2 plan row V3c).

`vector_step_jorba_zou(jets, eps_abs; eps_rel, vnorm)` is the norm-based
generalisation of `StepControl.step_jorba_zou` to a vector ODE: the
`k`-th Taylor coefficient is now a `d`-vector `c_k`, and the scalar
magnitude `|c_k|` is replaced by `vnorm(c_k)` (default the 2-norm).  The
opt-in adaptive policy `step_policy = :jorba_zou` wires this selector
into the `vector_solve_pade` driver.

These tests (`VSC.*`) are written RED-first per CLAUDE.md Rule 4; each
asserts an invariant against a known-correct value (Rule 5):

    VSC.1.1  d=1 reduction       — a 1-jet call equals StepControl
                                   .step_jorba_zou exactly
    VSC.1.2  known vector jet    — hand-computed coefficient norms ⇒
                                   closed-form min over k of (ε/‖c_k‖)^(1/k)
    VSC.1.3  norm choice         — 2-norm vs ∞-norm documented ordering
    VSC.1.4  adaptive driver     — :jorba_zou integrates the harmonic
                                   system, h-ceiling respected; :fixed
                                   output unchanged
    VSC.1.5  failure modes       — short jets, all-zero, ragged jets throw
    VSC.1.6  mutation-proof      — recorded at end of file

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/vector_step_control_test.jl`) and under
`runtests.jl`.
"""

using Test
using LinearAlgebra: norm
using PadeTaylor
using PadeTaylor.StepControl:       step_jorba_zou
using PadeTaylor.VectorStepControl: vector_step_jorba_zou
using PadeTaylor.VectorProblems:    VectorPadeTaylorProblem, vector_solve_pade

@testset "VectorStepControl (V3c)" begin

    @testset "VSC.1.1 d=1 reduces to StepControl.step_jorba_zou" begin
        # A single-jet call must reproduce the scalar selector EXACTLY
        # (the 2-norm of a length-1 vector is `abs` of its sole entry).
        # Several jets, including the canonical exp jet [1/k! for k=0..30].
        exp_jet = [1.0 / factorial(big(k)) for k in 0:30]
        exp_jet = Float64.(exp_jet)
        jets_to_test = [
            exp_jet,
            Float64[1.0, 0.5, 0.25, 0.125, 0.0625],
            Float64[2.0, -1.0, 0.3, -0.04, 0.005, -0.0006],
            Float64[1.0, 0.0, 0.0, 0.0, 7.0],   # exercises zero-skip
        ]
        for jet in jets_to_test
            for ε in (1e-12, 1e-8, 1e-6)
                scalar = step_jorba_zou(jet, ε)
                vector = vector_step_jorba_zou([jet], ε)
                @test vector == scalar
            end
        end
        # eps_rel passthrough: a non-default eps_rel must reduce too.
        jet = Float64[3.0, 1.0, 0.2, 0.01]
        @test vector_step_jorba_zou([jet], 1e-10; eps_rel = 1e-9) ==
              step_jorba_zou(jet, 1e-10; eps_rel = 1e-9)
    end

    @testset "VSC.1.2 known vector jet — closed-form norm formula" begin
        # Two components, order 3.  k-th coefficient vector
        #   c_k = [jet1[k+1], jet2[k+1]].
        # jet1 = [1, 1, 3, 4], jet2 = [0, 1, 4, 3].
        #   ‖c_2‖₂ = ‖[3,4]‖ = 5,  ‖c_3‖₂ = ‖[4,3]‖ = 5,  ‖c_0‖₂ = 1.
        # With p = 3, the formula uses k ∈ {2,3}:
        #   h = min((ε/5)^(1/2), (ε/5)^(1/3)).
        jet1 = Float64[1.0, 1.0, 3.0, 4.0]
        jet2 = Float64[0.0, 1.0, 4.0, 3.0]
        ε = 1e-9
        h = vector_step_jorba_zou([jet1, jet2], ε)
        expected = min((ε / 5.0)^(1 / 2), (ε / 5.0)^(1 / 3))
        @test h ≈ expected rtol = 1e-14
        @test h isa Float64

        # ε resolution: eps_rel branch.  ‖c_0‖₂ = ‖[6,8]‖ = 10, so with
        # eps_abs = 1e-12 < eps_rel·10 = 1e-9, ε = 1e-9 is used.
        j1 = Float64[6.0, 0.0, 3.0, 4.0]
        j2 = Float64[8.0, 0.0, 4.0, 3.0]
        h2 = vector_step_jorba_zou([j1, j2], 1e-12; eps_rel = 1e-10)
        εr = 1e-10 * 10.0
        @test h2 ≈ min((εr / 5.0)^(1 / 2), (εr / 5.0)^(1 / 3)) rtol = 1e-14
    end

    @testset "VSC.1.3 norm choice — 2-norm vs ∞-norm ordering" begin
        # For any vector, ‖v‖_∞ ≤ ‖v‖₂.  The selector divides ε by the
        # norm of the leading coefficient, so a SMALLER norm gives a
        # LARGER step.  Hence h(∞-norm) ≥ h(2-norm) on a jet whose
        # leading coefficients have ∞-norm strictly below their 2-norm.
        # jet1 = [1,0,3,3], jet2 = [0,0,3,3]: c_2 = c_3 = [3,3],
        #   ‖[3,3]‖₂ = 3√2 ≈ 4.243,  ‖[3,3]‖_∞ = 3.
        jet1 = Float64[1.0, 0.0, 3.0, 3.0]
        jet2 = Float64[0.0, 0.0, 3.0, 3.0]
        ε = 1e-9
        h_2   = vector_step_jorba_zou([jet1, jet2], ε)               # default 2-norm
        h_inf = vector_step_jorba_zou([jet1, jet2], ε;
                                      vnorm = v -> norm(v, Inf))
        @test h_inf > h_2
        # Closed-form check of both.
        @test h_2   ≈ min((ε / (3 * sqrt(2)))^(1 / 2),
                          (ε / (3 * sqrt(2)))^(1 / 3)) rtol = 1e-14
        @test h_inf ≈ min((ε / 3.0)^(1 / 2), (ε / 3.0)^(1 / 3)) rtol = 1e-14
    end

    @testset "VSC.1.4 adaptive driver — :jorba_zou vs :fixed" begin
        # Harmonic system y₁'=y₂, y₂'=−y₁, y0=[1,0] ⇒ y=[cos z,−sin z].
        harm = (z, y) -> [y[2], -y[1]]
        order = 24

        # :fixed default — byte-identical to the V3b behaviour.
        prob_f = VectorPadeTaylorProblem(harm, [1.0, 0.0], (0.0, 2.0); order)
        sol_fixed_default = vector_solve_pade(prob_f; h = 0.25)
        sol_fixed_explicit = vector_solve_pade(prob_f; h = 0.25,
                                               step_policy = :fixed)
        @test sol_fixed_default.z == sol_fixed_explicit.z
        @test sol_fixed_default.y == sol_fixed_explicit.y
        @test sol_fixed_default.h == sol_fixed_explicit.h

        # :jorba_zou — adaptive.  h = 0.25 is the CEILING; selected steps
        # must never exceed it (modulo the final clamp to z_end).
        prob_a = VectorPadeTaylorProblem(harm, [1.0, 0.0], (0.0, 2.0); order)
        sol_adaptive = vector_solve_pade(prob_a; h = 0.25,
                                         step_policy = :jorba_zou)
        for hk in sol_adaptive.h
            @test hk ≤ 0.25 + 1e-14
            @test hk > 0
        end
        @test sol_adaptive.z[end] ≈ 2.0
        # Same accuracy as :fixed against the closed form.
        for k in 1:length(sol_adaptive.z)
            z = sol_adaptive.z[k]
            @test sol_adaptive.y[k][1] ≈ cos(z)  atol = 1e-9
            @test sol_adaptive.y[k][2] ≈ -sin(z) atol = 1e-9
        end
        # Dense callable still works under the adaptive policy.
        for z in (0.3, 1.1, 1.7)
            yz = sol_adaptive(z)
            @test yz[1] ≈ cos(z)  atol = 1e-9
            @test yz[2] ≈ -sin(z) atol = 1e-9
        end

        # An unknown policy symbol must throw (Rule 1).
        @test_throws ArgumentError vector_solve_pade(
            prob_f; h = 0.25, step_policy = :bogus)
    end

    @testset "VSC.1.5 failure modes throw informatively" begin
        # Jets shorter than 2 — no notion of "last two coefficients".
        @test_throws ArgumentError vector_step_jorba_zou(
            [Float64[1.0]], 1e-9)
        @test_throws ArgumentError vector_step_jorba_zou(
            [Float64[]], 1e-9)

        # All coefficients zero — no nonzero coefficient to estimate from.
        @test_throws ArgumentError vector_step_jorba_zou(
            [Float64[0.0, 0.0, 0.0], Float64[0.0, 0.0, 0.0]], 1e-9)

        # Ragged jets of unequal length — a malformed vector jet.
        @test_throws ArgumentError vector_step_jorba_zou(
            [Float64[1.0, 0.5, 0.25], Float64[1.0, 0.5]], 1e-9)

        # Empty jet collection — d = 0, nothing to integrate.
        @test_throws ArgumentError vector_step_jorba_zou(
            Vector{Float64}[], 1e-9)
    end

    # -------------------------------------------------------------------------
    # VSC.1.6 — mutation-proof.  Three mutations were applied to
    # `src/VectorStepControl.jl` / `src/VectorProblems.jl`, the suite
    # re-run, RED counts recorded, mutations reverted.  Documented here;
    # this block is not executable (impl is in its restored state).
    # Baseline: 69/69 GREEN.
    #
    #   M1 — ignore the other components: in `_coef_vector`, return
    #        `[jets[1][k+1]]` instead of the full cross-component slice
    #        `[jets[i][k+1] for i in 1:d]`.  The selector then keys off
    #        component 1's magnitude alone — 4 instead of ‖[3,4]‖=5 in
    #        VSC.1.2, 0 in VSC.1.3's leading coefficients.  Observed:
    #        3 of 69 failed (1 in VSC.1.2, 2 in VSC.1.3) — BIT.
    #   M2 — wrong exponent: `(ε/‖c_k‖)^(1/k)` → `(ε/‖c_k‖)^(1/(k+1))`.
    #        Observed: 17 of 69 failed — every VSC.1.1 d=1-reduction
    #        assertion (13) plus VSC.1.2 (2) and VSC.1.3 (2).  The
    #        hardest-biting mutation: the exponent is load-bearing in
    #        every closed-form and reduction check — BIT.
    #   M3 — drop the h-ceiling cap in the :jorba_zou driver branch
    #        (`src/VectorProblems.jl`): use the raw `vector_step_jorba_zou`
    #        step `h_jz` uncapped instead of `min(h_jz, h_ceiling)`.  The
    #        harmonic jet's Jorba-Zou step at ε = eps(Float64) far
    #        exceeds the 0.25 ceiling, so the trajectory takes one giant
    #        step.  Observed: 6 of VSC.1.4's assertions failed (the
    #        `hk ≤ 0.25 + 1e-14` ceiling checks) — BIT.
    #
    # All three mutations bit; impl restored to the correct state and
    # both this suite (69/69) and the V3b suite (72/72) re-confirmed
    # GREEN.
    # -------------------------------------------------------------------------

end
