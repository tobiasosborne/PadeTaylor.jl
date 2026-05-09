"""
Phase 3 RED tests for `PadeTaylor.Coefficients.taylor_coefficients_1st`
and `taylor_coefficients_2nd` (per DESIGN.md §4 Phase 3).

Oracle source:
    external/probes/coefficients-oracle/capture.wl   (wolframscript primary)
    external/probes/coefficients-oracle/capture.py   (sympy/mpmath cross-check)
    external/probes/coefficients-oracle/verify.jl    (cross-validates and emits pinned file)
Pinned values: test/_oracle_coefficients.jl (auto-generated; do NOT hand-edit).

Test plan:
    3.1.1  f(z,y)=y, y(0)=1, order 10 Float64           — exp(z), spec-from-scratch
    3.1.2  f(z,y)=z²+y², y(0)=0, order 14 Float64       — Bessel-ratio oracle
    3.1.3  f(z,u,u')=6u², FW 2011 ICs, order 30 Float64 — Weierstrass oracle
    3.1.4  3.1.1 at order 60, BigFloat(precision=256)   — precision tier
    3.1.5  Mutation-proof procedure documented at end (executed manually after GREEN)

All Float64 assertions use a tight relative tolerance; the BigFloat-256
assertion uses 2⁻²⁰⁰ (matches the empirically-measured radius from
RESEARCH.md §3.3 probe, max 2.16e-78).
"""

using Test
using PadeTaylor.Coefficients: taylor_coefficients_1st, taylor_coefficients_2nd

include("_oracle_coefficients.jl")

@testset "Coefficients (Phase 3)" begin

    @testset "3.1.1 exp(z) Taylor, order 10, Float64" begin
        # f(z, y) = y, y(0) = 1 ⇒ y(z) = exp(z); c_k = 1/k!.
        # Oracle: wolframscript Series[Exp[z], {z, 0, 10}].
        # Cross-check: sympy Rational(1, factorial(k)).
        c = taylor_coefficients_1st((z, y) -> y, 0.0, 1.0, 10)
        @test length(c) == 11
        @test eltype(c) == Float64
        for k in 0:10
            @test isapprox(c[k+1], c_3_1_1_exp_o10[k+1]; rtol = 1e-14, atol = 1e-300)
        end
    end

    @testset "3.1.2 dy/dz = z² + y², y(0) = 0, order 14, Float64" begin
        # FW 2011 sec. 2.2.1 demo problem. Closed-form solution:
        #   y(t) = t · J_{3/4}(t²/2) / J_{-1/4}(t²/2)   (analytic at t=0).
        # Hand-verifiable leading coefs: y[3] = 1/3, y[7] = 1/63, y[11] = 2/2079.
        # Oracle: wolframscript Series of the closed form.
        # Cross-check: sympy explicit J_ν power-series expansion of the ratio.
        c = taylor_coefficients_1st((z, y) -> z^2 + y^2, 0.0, 0.0, 14)
        @test length(c) == 15
        @test eltype(c) == Float64
        for k in 0:14
            @test isapprox(c[k+1], c_3_1_2_bessel_o14[k+1]; rtol = 1e-13, atol = 1e-15)
        end
    end

    @testset "3.1.3 u'' = 6u², FW 2011 ICs, order 30, Float64" begin
        # Equianharmonic Weierstrass companion to PI. ICs string-matched to
        # FW 2011 markdown lines 292-295:
        #   c_1 = -1, c_2 = g_3 = 2,
        #   u(0)  = 1.071822516416917
        #   u'(0) = 1.710337353176786.
        # Oracle: wolframscript AsymptoticDSolveValue on the ODE with these ICs.
        # Cross-check: Series[WeierstrassP[z + c_1, {0, c_2}], {z, 0, 30}]
        # — the two methods agree to 3.78e-15 (verify.jl).
        u0  = 1.071822516416917
        up0 = 1.710337353176786
        c = taylor_coefficients_2nd((z, u, up) -> 6 * u^2, 0.0, u0, up0, 30)
        @test length(c) == 31
        @test eltype(c) == Float64
        # Coefficients grow ~linearly in k (the ODE is heading toward a pole),
        # so use relative tolerance throughout.
        for k in 0:30
            @test isapprox(c[k+1], c_3_1_3_weierstrass_o30[k+1]; rtol = 1e-12, atol = 1e-14)
        end
    end

    @testset "3.1.4 exp(z) Taylor, order 60, BigFloat(precision=256)" begin
        # Same ODE as 3.1.1, raised precision tier. Per RESEARCH.md §3.3
        # empirical probe, Taylor1 over a 256-bit substrate delivers > 200
        # genuine bits at order 60+ (probe measured max radius 2.16e-78 ≈
        # 2⁻²⁵⁸, well below the 2⁻²⁰⁰ tolerance asserted here).
        # Oracle: wolframscript N[1/k!, 80] for k=0..60.
        # Cross-check: mpmath at 80 decimal digits.
        setprecision(BigFloat, 256) do
            c = taylor_coefficients_1st((z, y) -> y,
                                         BigFloat(0), BigFloat(1), 60)
            @test length(c) == 61
            @test eltype(c) == BigFloat
            tol = BigFloat(2)^(-200)
            for k in 0:60
                expected = parse(BigFloat, c_3_1_4_exp_o60_str[k+1])
                Δ = abs(c[k+1] - expected)
                scale = max(abs(expected), BigFloat(2)^(-256))
                @test Δ / scale < tol
            end
        end
    end
end

#=
Phase 3 mutation-proof procedure (test 3.1.5) — VERIFIED 2026-05-09
──────────────────────────────────────────────────────────────────
Both the 1st- and 2nd-order recurrences were each mutation-tested
independently before commit; both mutations RED'd the test suite.

  Mutation A — off-by-one in the 1st-order integration step
  (`src/Coefficients.jl`, `taylor_coefficients_1st`):
      y[k] = f_t[k-1] / k     →     y[k] = f_t[k-1] / (k + 1)
    VERIFIED to fail 3.1.1 immediately (c[1] = 1.0 expected, ~0.5 got).
    Confirms the 1st-order recurrence is genuinely tested.

  Mutation C — off-by-one in the 2nd-order recurrence
  (`src/Coefficients.jl`, `taylor_coefficients_2nd`):
      u[j] = f_t[j-2] / (j*(j-1))   →   u[j] = f_t[j-2] / (j*(j+1))
    VERIFIED to fail 3.1.3 at c[2] onward.
    Confirms the 2nd-order recurrence is genuinely tested (Mutation A
    by itself only exercises the 1st-order code path; 3.1.4 is also
    1st-order at higher precision, so without Mutation C the entire
    `taylor_coefficients_2nd` body would be unproven).

Per worklog 001 ("Frictions surfaced — mutation-proof procedure was
wrong on first attempt") and CLAUDE.md Rule 4 + Rule 5: each load-bearing
recurrence got an independent mutation that REDs a test before this
commit.
=#
