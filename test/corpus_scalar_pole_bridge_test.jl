# test/corpus_scalar_pole_bridge_test.jl
#
# ============================================================================
# Ground-truth corpus: SCALAR pole-bridging (bead padetaylor-2tj5,
# epic padetaylor-p1v0).
#
# WHAT THIS FILE PINS
# -------------------
# The whole point of PadeTaylor.jl is in README §2: a ratio of polynomials
# can see *past* a pole where plain Taylor truncation cannot.  This file
# nails that claim down against eight independent closed-form / mpmath
# oracles, each a scalar ODE whose exact solution is known across a pole:
#
#   CPB.1  wp-equianharmonic    u''=6u^2,  u=WP(z-1;0,2)   double pole z=1
#   CPB.2  tan (1st->2nd order) u''=2u+2u^3, u=tan(z)      simple pole pi/2
#   CPB.3  simple-pole-riccati  u''=2u^3,    u=1/(1-z)     simple pole z=1
#   CPB.4  coth-real-pole       u''=-2u+2u^3, u=coth(z)    simple pole z=0
#   CPB.5  fw-riccati-bessel    u''=2t+2t^2 u+2u^3         pole t~2.0031
#   CPB.6  tanh-imaginary-pole  u''=-2u+2u^3, u=tanh(z)    pole i*pi/2
#   CPB.7  jacobi-sn complex    u''=-(5/4)u+(1/2)u^3       pole i*Kp
#   CPB.8  Taylor-foil keystone Padé bridges; Taylor diverges  (the headline)
#
# GROUND TRUTH each case embodies
# -------------------------------
# All recipes are verbatim from docs/test_corpus/01_corpus_catalogue.md
# "Full candidate records".  The first-order Riccati ODEs (tan, simple-pole,
# coth, tanh, fw) are recast as SECOND-order ODEs by one differentiation,
# because solve_pade's public API only accepts the 2nd-order (tuple-IC)
# branch (src/Problems.jl:178); the recast is exact:
#   u'=1+u^2  => u''=2u u'=2u(1+u^2)=2u+2u^3           (tan)
#   u'=u^2    => u''=2u u'=2u^3                         (simple pole)
#   u'=1-u^2  => u''=-2u u'=-2u(1-u^2)=-2u+2u^3         (coth, tanh)
#   u'=t^2+u^2=> u''=2t+2u u'=2t+2t^2 u+2u^3            (fw-riccati-bessel)
# The recasts are themselves catalogue records (tan-2nd-order, etc).
#
# Real-axis poles use solve_pade with a single segment that brackets the
# pole strictly inside [0,1] in the rescaled variable (the worklog-004
# idiom reproduced in problems_test.jl).  Off-real-axis poles (tanh, sn)
# use path_network_solve + eval_at, the sanctioned complex-direction path
# (src/PathNetwork.jl), since solve_pade's `min(h_max, z_end-z)` is real-only.
#
# TOLERANCES are justified by the METHOD's accuracy, never fitted to pass:
#   * Float64 (15,15)-diagonal Padé past a pole saturates at a ~3e-10
#     rational-approximation floor at t~0.7 (README §2; problems_test.jl
#     6.1.3 uses rtol=1e-9; 6.1.6 docstring quantifies the 3.5e-10 floor).
#     We use rtol 1e-9 for first-order-pole bridges that land near t~0.7.
#   * Exact-rational solutions (1/(1-z)) hit machine precision: rtol 1e-13.
#   * Complex single-step bridges with large |h| (sn h=Kp~2.16, tanh h=pi)
#     accumulate more rational-approx error: 1e-7 (tanh), 1e-5 (sn, h=Kp
#     spans a full half-period).  Documented per-assertion below.
#   * fw past a NON-algebraic Bessel pole in one h=2.5 step: 1e-5.
#
# REFERENCES: FW2011 §2.2.1 (fw-riccati-bessel, md:122-134) and §5.1
# (WP keystone, md:281-302); DLMF 4.21 (tan), 4.28 (coth/tanh),
# 22.13.1 (sn).  Oracle constants captured by
# external/probes/corpus-oracles/scalar-pole-bridge/capture.py at dps=50.
#
# MUTATION-PROOF (Rule 4): see the comment block at the end of the file —
# each load-bearing assertion was verified to go RED under a deliberate
# oracle perturbation, then restored.
# ============================================================================

using Test
using PadeTaylor
import PadeTaylor: taylor_coefficients_2nd

include(joinpath(@__DIR__, "_oracle_corpus_scalar_pole_bridge.jl"))

@testset "Corpus: scalar pole-bridge (CPB)" begin

    # ----------------------------------------------------------------------
    # CPB.1  wp-equianharmonic : u''=6u^2, double pole at z=1.
    # ----------------------------------------------------------------------
    @testset "CPB.1  WP(z-1;0,2): Padé bridges the double pole at z=1" begin
        fW(z, u, up) = 6 * u^2
        prob = PadeTaylorProblem(fW, (WP_U0, WP_UP0), (0.0, 1.5); order = 30)
        sol  = solve_pade(prob; h_max = 1.5)         # single seg, pole at t~0.667
        u_before, _ = sol(0.5)
        u_past, up_past = sol(1.05)                  # PAST the pole
        u_far, _ = sol(1.4)
        # rtol justified above: Float64 (15,15) Padé past-pole floor ~3e-10.
        @test isapprox(u_before, WP_U_AT_0p5;  rtol = 1e-12)
        @test isapprox(u_past,   WP_U_AT_1p05; rtol = 1e-9)
        @test isapprox(up_past,  WP_UP_AT_1p05; rtol = 1e-9)
        @test isapprox(u_far,    WP_U_AT_1p4;  rtol = 1e-7)
    end

    @testset "CPB.1.bf  WP at BigFloat-256: tighter past-pole accuracy" begin
        setprecision(BigFloat, 256) do
            fW(z, u, up) = 6 * u^2
            u0  = parse(BigFloat, "1.07182251641691742068161429481250292922453268")
            up0 = parse(BigFloat, "1.71033735317678620846433189149624950839231872")
            zsp = (BigFloat(0), BigFloat(3) / 2)
            prob = PadeTaylorProblem(fW, (u0, up0), zsp; order = 40)
            sol  = solve_pade(prob; h_max = BigFloat(3) / 2)
            u_ref = parse(BigFloat, WP_U_AT_1p05_STR)
            u_eval, _ = sol(BigFloat(105) / 100)
            # order=40 BF-256 converges to the IC-precision floor (see
            # problems_test.jl 6.1.6); rtol 1e-13 is well above it.
            @test isapprox(u_eval, u_ref; rtol = BigFloat("1e-13"))
        end
    end

    # ----------------------------------------------------------------------
    # CPB.2  tan : u''=2u+2u^3, simple pole at pi/2, bridge h=3 (pi/2 interior).
    # ----------------------------------------------------------------------
    @testset "CPB.2  tan(z): bridge the simple pole at pi/2" begin
        ftan(z, u, up) = 2 * u + 2 * u^3
        prob = PadeTaylorProblem(ftan, (0.0, 1.0), (0.0, 3.0); order = 30)
        sol  = solve_pade(prob; h_max = 3.0)         # pi/2 at t~0.524, interior
        u2, up2 = sol(2.0)                           # PAST pole pi/2
        u25, _  = sol(2.5)
        @test isapprox(u2,  TAN_U_AT_2p0;  rtol = 1e-9)
        @test isapprox(up2, TAN_UP_AT_2p0; rtol = 1e-9)  # u' = 1+tan^2 = sec^2
        @test isapprox(u25, TAN_U_AT_2p5;  rtol = 1e-9)
    end

    # ----------------------------------------------------------------------
    # CPB.3  simple-pole-riccati : u''=2u^3, u=1/(1-z), pole z=1 exact rational.
    # ----------------------------------------------------------------------
    @testset "CPB.3  1/(1-z): exact-rational bridge of the pole at z=1" begin
        fsp(z, u, up) = 2 * u^3
        prob = PadeTaylorProblem(fsp, (1.0, 1.0), (0.0, 2.0); order = 30)
        sol  = solve_pade(prob; h_max = 2.0)         # pole z=1 at t=0.5 (midpoint)
        u_before, _ = sol(0.5)
        u_past, up_past = sol(1.5)                   # PAST pole
        u_end, _ = sol(2.0)
        # Exact rationals: u(0.5)=2, u(1.5)=-2, u(2)=-1; up = u^2.
        @test isapprox(u_before, 2.0;  rtol = 1e-13)
        @test isapprox(u_past,  -2.0;  rtol = 1e-13)
        @test isapprox(up_past,  4.0;  rtol = 1e-12)  # u' = u^2 = (-2)^2 = 4
        @test isapprox(u_end,   -1.0;  rtol = 1e-13)
    end

    # ----------------------------------------------------------------------
    # CPB.4  coth-real-pole : u''=-2u+2u^3, simple pole z=0, bridge [-1,1].
    # ----------------------------------------------------------------------
    @testset "CPB.4  coth(z): bridge the real-axis pole at z=0" begin
        fcoth(z, u, up) = -2 * u + 2 * u^3
        prob = PadeTaylorProblem(fcoth, (COTH_U0_M1, COTH_UP0_M1),
                                 (-1.0, 1.0); order = 30)
        sol  = solve_pade(prob; h_max = 2.0)         # pole z=0 at t=0.5 (midpoint)
        u1, up1 = sol(1.0)                           # PAST pole; odd symm u(1)=-u(-1)
        @test isapprox(u1,  COTH_U_AT_1p0;  rtol = 1e-9)
        @test isapprox(up1, COTH_UP_AT_1p0; rtol = 1e-8)
        # Odd-symmetry invariant: coth(1) = -coth(-1).
        @test isapprox(u1, -COTH_U0_M1; rtol = 1e-9)
    end

    # ----------------------------------------------------------------------
    # CPB.5  fw-riccati-bessel : the FW2011 §2.2.1 motivating demo.  Pole at
    # the NON-algebraic Bessel zero t1=2.0031; Padé steps straight through it.
    # ----------------------------------------------------------------------
    @testset "CPB.5  FW Riccati-Bessel: bridge the t~2.003 Bessel pole" begin
        ffw(z, u, up) = 2 * z + 2 * z^2 * u + 2 * u^3
        prob = PadeTaylorProblem(ffw, (0.0, 0.0), (0.0, 2.5); order = 30)
        sol  = solve_pade(prob; h_max = 2.5)         # t1~2.003 at t~0.801, interior
        u05, _ = sol(0.5)
        u15, _ = sol(1.5)                            # just BEFORE the pole
        u21, _ = sol(2.1)                            # PAST the pole
        @test isapprox(u05, FW_U_AT_0p5; rtol = 1e-10)
        @test isapprox(u15, FW_U_AT_1p5; rtol = 1e-9)
        # Past a NON-algebraic pole in one big step: 1e-5 (justified above).
        @test isapprox(u21, FW_U_AT_2p1; rtol = 1e-5)
        # The pole really is bridged: u(2.1) is finite and the right sign.
        @test isfinite(u21) && u21 < 0
    end

    # ----------------------------------------------------------------------
    # CPB.6  tanh-imaginary-pole : complex-direction bridge through i*pi/2.
    # ----------------------------------------------------------------------
    @testset "CPB.6  tanh(z): imaginary-axis bridge through pole at i*pi/2" begin
        ftanh(z, u, up) = -2 * u + 2 * u^3
        prob = PadeTaylorProblem(ftanh,
                                 (complex(0.0, 0.0), complex(1.0, 0.0)),
                                 (complex(0.0, 0.0), complex(0.0, 1.0)); order = 30)
        ztarget = complex(0.0, 3 * pi / 4)           # PAST pole at i*pi/2
        sol = path_network_solve(prob, [ztarget]; h = float(pi), order = 30)
        u_past, up_past = eval_at(sol, ztarget)
        u_pre, _ = eval_at(sol, complex(0.0, pi / 4))   # BEFORE pole: tanh=i exact
        # tanh(i*3pi/4) = -i exact, tanh(i*pi/4) = i exact (catalogue).
        @test isapprox(u_past, complex(0.0, -1.0); atol = 1e-7)
        @test isapprox(up_past, complex(2.0, 0.0); atol = 1e-7)  # u'=1-u^2=1-(-i)^2=2
        @test isapprox(u_pre,  complex(0.0,  1.0); atol = 1e-9)
    end

    # ----------------------------------------------------------------------
    # CPB.7  jacobi-sn complex bridge : exact algebraic ICs, pole at i*Kp.
    # ----------------------------------------------------------------------
    @testset "CPB.7  sn(z,1/4): complex bridge through the pole at i*Kp" begin
        fsn(z, u, up) = -(5 // 4) * u + (1 // 2) * u^3
        uic  = complex(0.0, sqrt(2.0))               # u(i*Kp/2) = i*sqrt(2) exact
        upic = complex(3 / sqrt(2.0), 0.0)           # u'(i*Kp/2) = 3/sqrt(2) exact
        zic  = complex(0.0, SN_KP / 2)
        ztgt = complex(0.0, 3 * SN_KP / 2)           # PAST pole at i*Kp
        prob = PadeTaylorProblem(fsn, (uic, upic), (zic, ztgt); order = 30)
        sol  = path_network_solve(prob, [ztgt]; h = SN_KP, order = 30)
        u_past, _ = eval_at(sol, ztgt)
        # u(i*3Kp/2) = -i*sqrt(2) exact.  rtol 1e-5: one h=Kp (~2.16) step
        # spans a full imaginary half-period across the pole.
        @test isapprox(u_past, complex(0.0, -sqrt(2.0)); rtol = 1e-5)
    end

    # ----------------------------------------------------------------------
    # CPB.8  THE HEADLINE: same Taylor coefficients, Padé bridges, Taylor
    # diverges.  README §2 / problems_test.jl 6.1.5 generalised to tan + fw.
    # ----------------------------------------------------------------------
    @testset "CPB.8  Padé bridges where plain Taylor diverges" begin
        # tan: pole at pi/2; evaluate at z=2.0, past the natural radius.
        ftan(z, u, up) = 2 * u + 2 * u^3
        prob = PadeTaylorProblem(ftan, (0.0, 1.0), (0.0, 3.0); order = 30)
        sol  = solve_pade(prob; h_max = 3.0)
        u_pade, _ = sol(2.0)
        c_tan = taylor_coefficients_2nd(ftan, 0.0, 0.0, 1.0, 30)
        u_taylor = PadeTaylor.taylor_eval(c_tan, 2.0)
        @test isapprox(u_pade, TAN_U_AT_2p0; rtol = 1e-9)            # Padé wins
        @test abs(u_taylor - TAN_U_AT_2p0) / abs(TAN_U_AT_2p0) > 1.0 # Taylor loses

        # fw-riccati-bessel: pole at t~2.003; evaluate at z=2.1.
        ffw(z, u, up) = 2 * z + 2 * z^2 * u + 2 * u^3
        probf = PadeTaylorProblem(ffw, (0.0, 0.0), (0.0, 2.5); order = 30)
        solf  = solve_pade(probf; h_max = 2.5)
        uf_pade, _ = solf(2.1)
        c_fw = taylor_coefficients_2nd(ffw, 0.0, 0.0, 0.0, 30)
        uf_taylor = PadeTaylor.taylor_eval(c_fw, 2.1)
        @test isapprox(uf_pade, FW_U_AT_2p1; rtol = 1e-5)            # Padé wins
        @test abs(uf_taylor - FW_U_AT_2p1) / abs(FW_U_AT_2p1) > 0.1  # Taylor loses
    end
end

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
#
# M-oracle (CPB.3): in this file, assert `u_past == +2.0` instead of `-2.0`
#   (the exact-rational value 1/(1-1.5) = -2).  Result: CPB.3 went RED at
#   the flipped line (got -2.0, expected +2.0) and ONLY there; the other 26
#   assertions stayed GREEN.  Restored the `-2.0` literal; suite GREEN.
#
# M-impl (keystone, the gold standard): in src/PadeStepper.jl `_evaluate_pade`,
#   multiply the numerator constant term by 1.0000001 (a 1e-7 relative
#   perturbation of the stored Padé).  Result: 8 value-assertions went RED —
#   CPB.1 (u, u_past), CPB.1.bf, CPB.3 (3 assertions), CPB.4 (u, odd-symm) —
#   i.e. every assertion routed through the stored Padé caught the regression.
#   Restored src/PadeStepper.jl byte-for-byte (verified via `diff`); suite GREEN.
#
# COVERAGE NOTE worth recording: an earlier M-impl that perturbed the
# CHAIN-RULE factor in `pade_step_with_pade!` (`new_up / h_T -> / h_T^2`) did
# NOT redden any CPB assertion.  Reason: `sol(z)` derivatives are computed by
# `Problems.PadeTaylorSolution`'s OWN chain rule on the stored Padé
# (src/Problems.jl:235), not by the stepper's internal `state.up`.  So the
# `up`/derivative assertions here pin `Problems.jl`'s chain rule, while the
# stepper's `state.up` is exercised by the multi-segment propagation tests
# elsewhere (padestepper_test.jl).  This is a true coverage fact, not a gap
# in these single-segment cases.
#
# Run standalone with:
#   julia --project=. test/corpus_scalar_pole_bridge_test.jl
# ============================================================================
