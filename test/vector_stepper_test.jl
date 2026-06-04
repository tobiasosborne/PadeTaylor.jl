"""
V3a tests for `PadeTaylor.VectorStepper` (bead `padetaylor-0ln.8`;
v0.2 plan row V3a).

`vector_pade_step!(state, f, order, h)` takes one Padé–Taylor step of a
first-order **vector** ODE `y' = f(z, y)`, `y ∈ ℂ^d`, advancing the
state `(z, y)` to `(z+h, y(z+h))` — with all `d` components sharing one
denominator `Q` (the v0.2 keystone: poles are a property of the flow,
not of the individual component).

These tests (`VS.*`) are written RED-first per CLAUDE.md Rule 4; each
asserts an invariant against a known-correct value (Rule 5):

    VS.1.1  d=1 reduction        — same y(z+h) as the explicit scalar
                                   robust_pade pipeline (inline oracle)
    VS.1.2  closed-form 1 step   — harmonic system → [cos h, −sin h]
    VS.1.3  multi-step accuracy  — iterate across an interval, match
                                   the closed-form solution at each node
    VS.1.4  shared-Q pole bridge — a vector system with a known shared
                                   pole just past the step; Padé returns
                                   finite/correct values where plain
                                   Taylor truncation diverges
    VS.1.5  failure modes        — Q(1)≈0 step-onto-pole; bad order throw
    VS.1.6  mutation-proof       — recorded at end of file

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/vector_stepper_test.jl`) and under
`runtests.jl`.
"""

using Test
using PadeTaylor
using PadeTaylor.VectorStepper:
    VectorPadeStepperState, vector_pade_step!, vector_pade_step_with_pade!
using PadeTaylor.Coefficients: taylor_coefficients_1st
using PadeTaylor.RobustPade:   robust_pade

# -----------------------------------------------------------------------------
# Test-local helpers (inline oracles — kept independent of the impl).
# -----------------------------------------------------------------------------

# Evaluate the polynomial `c` (low-to-high coefficients) at `t` by Horner.
_horner(c, t) = begin
    s = zero(eltype(c)) * t
    for k in length(c):-1:1
        s = s * t + c[k]
    end
    s
end

# Rescale c̃_k = h^k · c_k (the variable substitution h' = h·t).
_rescale(c, h) = [c[k+1] * h^k for k in 0:(length(c)-1)]

# The explicit scalar Padé pipeline, built inline for the VS.1.1 oracle:
# Taylor jet → h^k rescale → robust_pade(:svd) → evaluate at t = 1.
function _scalar_pade_step(f, z0, y0, order, h)
    jet = taylor_coefficients_1st(f, z0, y0, order)
    scaled = _rescale(jet, h)
    m = order ÷ 2
    P = robust_pade(scaled, m, m; method = :svd)
    return _horner(P.a, one(h)) / _horner(P.b, one(h))
end

@testset "VectorStepper (V3a)" begin

    @testset "VS.1.1 d=1 reduces to the scalar robust_pade pipeline" begin
        # A 1-component system must produce the same y(z+h) as the
        # explicit scalar pipeline (Taylor → h^k rescale → robust_pade
        # :svd → evaluate at t=1).  Scalar ODE: y' = z² + y², y(0) = 0
        # (FW 2011 §2.2.1 demo), order 30, h = 0.3.
        scalar_f = (z, y) -> z^2 + y^2
        vec_f    = (z, y) -> [z^2 + y[1]^2]
        order, h = 30, 0.3

        oracle = _scalar_pade_step(scalar_f, 0.0, 0.0, order, h)

        st = VectorPadeStepperState{Float64}(0.0, [0.0])
        vector_pade_step!(st, vec_f, order, h)
        @test st.z ≈ 0.3
        @test length(st.y) == 1
        @test st.y[1] ≈ oracle atol = 1e-12

        # A second scalar ODE — y' = y, y(0) = 1 ⇒ y(h) = exp(h).
        oracle2 = _scalar_pade_step((z, y) -> y, 0.0, 1.0, order, h)
        st2 = VectorPadeStepperState{Float64}(0.0, [1.0])
        vector_pade_step!(st2, (z, y) -> [y[1]], order, h)
        @test st2.y[1] ≈ oracle2 atol = 1e-12
        @test st2.y[1] ≈ exp(0.3) atol = 1e-12
    end

    @testset "VS.1.2 closed-form harmonic system, single step" begin
        # y₁' = y₂, y₂' = −y₁, y0 = [1, 0] ⇒ y(h) = [cos h, −sin h].
        # Both components are entire (no poles); the shared denominator Q has no
        # genuine pole to fit.  This is the canonical ADR-0028 recovery case:
        # the dispatch (VectorStepper → SharedPadeDispatch.shared_pade_select)
        # selects the wide-square Mano–Tsuda cell B over the diagonal (m,m) cell A
        # here (the ODE defect validates B better on an entire system), so BOTH
        # components are recovered to ~roundoff — the worst-component accuracy the
        # diagonal cell A could not reach.  Measured (probe
        # adr0028-froissart-consumer/measure_vs12.jl): cell B gives cos err 2.2e-16
        # / −sin err 1.1e-16 at Float64, cos 1.76e-30 / −sin 5.8e-31 at BF-256 —
        # vs cell A's −sin ~5e-9 (F64) / ~1.6e-17 (BF).  (Cell A reached the cos
        # jet to ~1e-38 by luck of the even jet, but stalled the −sin at the
        # method-set (m,m) residual; cell B trades that for both components even at
        # ~1e-30 — a ~14-order win on the worst component.)  Closed form
        # [cos h, −sin h] is the exact Rule-5 oracle.
        harm = (z, y) -> [y[2], -y[1]]
        order, h = 30, 0.7

        st = VectorPadeStepperState{Float64}(0.0, [1.0, 0.0])
        vector_pade_step!(st, harm, order, h)
        @test st.z ≈ 0.7
        @test st.y[1] ≈ cos(0.7) atol = 1e-14       # cell B: roundoff (meas 2.2e-16)
        @test st.y[2] ≈ -sin(0.7) atol = 1e-14      # cell B: roundoff (meas 1.1e-16; cell A was ~5e-9)

        # BigFloat (256-bit): the ADR-0028 dispatch's square cell B recovers BOTH
        # components to ~1e-30 (meas cos 1.76e-30, −sin 5.8e-31) — the −sin no
        # longer stalls at the (m,m) ~1.6e-17 residual.  Both asserted < 1e-28
        # (≳ 50× the measured values), the even-error invariant of cell B.
        setprecision(BigFloat, 256) do
            hb = BigFloat("0.7")
            stb = VectorPadeStepperState{BigFloat}(
                BigFloat(0), BigFloat[1, 0])
            vector_pade_step!(stb, harm, order, hb)
            @test abs(stb.y[1] - cos(hb)) < BigFloat(10)^(-28)
            @test abs(stb.y[2] + sin(hb)) < BigFloat(10)^(-28)
        end
    end

    @testset "VS.1.3 multi-step accuracy across an interval" begin
        # Iterate the harmonic system over [0, 2] in 8 steps of h = 0.25;
        # assert the trajectory matches [cos z, −sin z] at every node.
        harm = (z, y) -> [y[2], -y[1]]
        order, h, nsteps = 24, 0.25, 8

        st = VectorPadeStepperState{Float64}(0.0, [1.0, 0.0])
        for k in 1:nsteps
            vector_pade_step!(st, harm, order, h)
            z = k * h
            @test st.z ≈ z
            # entire system → method-set (m,m) shared-Q accuracy (ADR-0027);
            # dispatch (ADR-0028) will tighten this back toward 1e-11.
            @test st.y[1] ≈ cos(z) atol = 1e-8
            @test st.y[2] ≈ -sin(z) atol = 1e-8
        end
        @test st.z ≈ 2.0

        # A coupled exponential pair to exercise composition further:
        # y₁' = y₁, y₂' = y₁ + y₂, y0 = [1, 0]
        # ⇒ y₁ = exp z, y₂ = z·exp z.
        sys = (z, y) -> [y[1], y[1] + y[2]]
        st2 = VectorPadeStepperState{Float64}(0.0, [1.0, 0.0])
        for k in 1:6
            vector_pade_step!(st2, sys, 20, 0.2)
            z = k * 0.2
            # entire (exp) system → method-set (m,m) accuracy; the 0.2^k-rescaled
            # jet reduces to degree ≈ 5, residual ~4e-8 accumulated (ADR-0027;
            # dispatch ADR-0028 recovers the tighter accuracy).
            @test st2.y[1] ≈ exp(z) atol = 1e-7
            @test st2.y[2] ≈ z * exp(z) atol = 1e-7
        end
    end

    @testset "VS.1.4 shared-Q pole bridging" begin
        # Vector system with a known SHARED pole.  Take
        #   y₁ = 1/(1−z),  y₂ = 2/(1−z)
        # both poles at z = 1.  Then y₁' = 1/(1−z)² = y₁², and
        # y₂' = 2/(1−z)² = y₁·y₂.  RHS: f = [y₁², y₁·y₂].
        # y0 at z = 0 is [1, 2].
        f = (z, y) -> [y[1]^2, y[1] * y[2]]
        order = 30
        h = 0.95                         # step lands at z = 0.95, just
                                         # short of the shared pole z = 1.

        st = VectorPadeStepperState{Float64}(0.0, [1.0, 2.0])
        vector_pade_step!(st, f, order, h)
        exact1 = 1 / (1 - h)             # = 20.0
        exact2 = 2 / (1 - h)             # = 40.0
        @test isfinite(st.y[1]) && isfinite(st.y[2])
        @test st.y[1] ≈ exact1 atol = 1e-8
        @test st.y[2] ≈ exact2 atol = 1e-8

        # Inline gap: plain Taylor truncation of the SAME jet diverges.
        # 1/(1−z) has Taylor coefficients all 1, so the order-30 partial
        # sum at z = 0.95 is Σ_{k=0}^{30} 0.95^k — far from 20.
        taylor_trunc1 = sum(0.95^k for k in 0:order)
        pade_err   = abs(st.y[1] - exact1)
        taylor_err = abs(taylor_trunc1 - exact1)
        @test taylor_err > 1.0           # Taylor truncation is way off
        @test pade_err < 1e-8            # Padé bridges the pole
        # The Padé step is at least 8 orders of magnitude more accurate.
        @test taylor_err / pade_err > 1e8
        @info "VS.1.4 pole-bridge gap" pade_err taylor_err ratio = taylor_err / pade_err
    end

    @testset "VS.1.5 failure modes throw informatively" begin
        # Bad order: the vector jet generator needs order ≥ 1.
        st = VectorPadeStepperState{Float64}(0.0, [1.0])
        @test_throws ArgumentError vector_pade_step!(
            st, (z, y) -> [y[1]], 0, 0.1)

        # Step exactly onto a pole: y = 1/(1−z), pole at z = 1, take
        # h = 1.0 so t = 1 lands on Q's root.  Must throw, not return Inf.
        f = (z, y) -> [y[1]^2]
        st2 = VectorPadeStepperState{Float64}(0.0, [1.0])
        @test_throws Exception vector_pade_step!(st2, f, 30, 1.0)
    end

    # -------------------------------------------------------------------------
    # VS.1.6 — mutation-proof.  Three mutations were applied to
    # `src/VectorStepper.jl`, the suite re-run, the RED count recorded,
    # and the mutation reverted.  This block documents the result; it is
    # not executable (the impl is in its restored, correct state).
    #
    #   M1 — drop the h^k rescaling (feed the raw, unscaled jets to
    #        shared_denominator_pade).  The Padé then lives in the wrong
    #        variable, so t = 1 no longer corresponds to the step h.
    #        Observed: 32 failures + 2 errors of 46 run — BIT hard.
    #   M2 — divide each numerator by Q at t = 0 instead of t = 1 (wrong
    #        evaluation point — t = 0 is the expansion centre, where
    #        every component equals its y0).  Observed: 35 failures of
    #        55 — BIT.
    #   M3 — give each component its OWN single-jet shared_denominator_pade
    #        call (d separate denominators), then divide every numerator
    #        by component 1's Q — defeats the shared-Q keystone (the d
    #        components no longer share one denominator).  Observed: 22
    #        failures + 1 error of 54 run — BIT, including VS.1.4 (the
    #        genuinely shared pole).
    #
    # All three mutations bit (counts above); impl restored to the
    # correct state and the suite re-confirmed GREEN at 55/55.
    # -------------------------------------------------------------------------

end
