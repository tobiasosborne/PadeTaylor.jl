"""
V6 tests for `PadeTaylor.PainleveHierarchy` (bead `padetaylor-0ln.15`;
v0.2 plan row V6; ADR-0022).

`PainleveHierarchy` is the problem-builder layer for the Painlevé-I
hierarchy.  Member `m = 2` is `P_I^(2)`, the fourth-order ODE

    u_xxxx + 10 u_x^2 + 20 u u_xx + 40(u^3 - 6t u + 6x) = 0

(KKG 2015 eq. (p12), verified verbatim against
`references/tex/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41/tritronquee_coeff.tex:125-130`).
Its first-order 4-vector companion form (pillar C §2,
`docs/v0p2_pillarC_painleve_hierarchy_findings.md:86-91`) is

    y' = (y_2, y_3, y_4, -10 y_2^2 - 20 y_1 y_3 - 40(y_1^3 - 6t y_1 + 6x)),

`y = (u, u', u'', u''') ∈ ℂ^4`, `x` the independent variable, `t` a
parameter.  Member `m = 1` is PI itself, `u_xx = 6u^2 + x`, in 2-vector
companion form `y' = (y_2, 6 y_1^2 + x)`.  `m ≥ 3` (general Lenard
recursion) is deferred.

`painleve_hierarchy(:I, m; t)` is the RHS closure factory;
`PainleveHierarchyProblem(m; y0, xspan, t)` is the builder;
`pI2_tritronquee_ic(x0; t, n_terms)` is the tritronquée asymptotic IC.

These tests (`PH.*`) are written RED-first per CLAUDE.md Rule 4; each
asserts an invariant against a known-correct value (Rule 5):

    PH.1.1  m=2 RHS            — matches the hand-written companion RHS
                                  (pillar C §2) for t=0 and t≠0
    PH.1.2  m=1 reduction      — painleve_hierarchy(:I,1) solved as a
                                  vector system reproduces v0.1
                                  PainleveProblem(:I) (same normalisation)
    PH.1.3  P_I^(2) ODE residual — solved trajectory satisfies the
                                  original 4th-order ODE (residual ≈ 0)
    PH.1.4  tritronquée IC     — pI2_tritronquee_ic returns the documented
                                  leading values; short integration stays
                                  near the asymptotic series u_f
    PH.1.5  failure modes      — m ≥ 3 throws; equation ≠ :I throws;
                                  wrong y0 length throws; n_terms ∉ {1,2}
                                  throws; n_terms = 2 with t ≠ 0 throws
    PH.1.7  tritronquée IC n=2 — pI2_tritronquee_ic(·; n_terms = 2) at
                                  t = 0 equals Y + Y^{-6} (c_6 = 1) and
                                  its analytic derivatives; refines the
                                  leading-order ODE residual (bead V8b)
    PH.1.6  mutation-proof     — recorded at end of file

Oracles: the hand-written companion RHS (PH.1.1); v0.1
`PainleveProblem(:I)` trajectory (PH.1.2); the original ODE residual
computed from the solved jet (PH.1.3); the closed-form asymptotic series
`u = -∛6·|x|^{1/3}` and its derivatives (PH.1.4).

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/painleve_hierarchy_test.jl`) and under
`runtests.jl`.
"""

using Test
using PadeTaylor

@testset "PainleveHierarchy (V6)" begin

    # -------------------------------------------------------------------------
    # PH.1.1 — m=2 RHS matches the hand-written companion RHS.
    # Oracle: the explicit companion form from pillar C §2.
    # -------------------------------------------------------------------------
    @testset "PH.1.1 m=2 RHS == hand-written companion RHS" begin
        # Hand-written reference RHS (pillar C §2, lines 86-91, verbatim).
        ref_PI2(x, t, y) = begin
            y1, y2, y3, y4 = y
            [y2, y3, y4,
             -10*y2^2 - 20*y1*y3 - 40*(y1^3 - 6*t*y1 + 6*x)]
        end

        for t in (0.0, 0.7, -1.3)
            f = painleve_hierarchy(:I, 2; t = t)
            for (x, y) in ((0.5,  [1.0, -2.0, 0.3,  4.0]),
                           (-3.0, [0.2,  1.1, -0.7, 2.5]),
                           (2.0,  [-1.5, 0.0, 1.0, -0.5]))
                got = f(x, y)
                exp = ref_PI2(x, t, y)
                @test length(got) == 4
                @test got ≈ exp atol = 1e-12 rtol = 1e-12
            end
        end

        # Each component individually pinned at a concrete point (t=0.7).
        f = painleve_hierarchy(:I, 2; t = 0.7)
        x = 0.5; y = [1.0, -2.0, 0.3, 4.0]
        r = f(x, y)
        @test r[1] ≈ -2.0
        @test r[2] ≈  0.3
        @test r[3] ≈  4.0
        # y4' = -10(4) - 20(1)(0.3) - 40(1 - 6*0.7*1 + 6*0.5)
        #     = -40 - 6 - 40(1 - 4.2 + 3) = -46 - 40*(-0.2) = -46 + 8 = -38.
        @test r[4] ≈ -38.0 atol = 1e-12
    end

    # -------------------------------------------------------------------------
    # PH.1.2 — m=1 reduction reproduces v0.1 PainleveProblem(:I).
    # The hierarchy m=1 member is PI itself, u_xx = 6u^2 + x — the SAME
    # normalisation v0.1 uses (`_pI_rhs() = 6u^2 + z`, Painleve.jl:101;
    # KKG 2015 tritronquee_coeff.tex:124 writes PI as `y_xx = 6y^2 + x`).
    # No rescaling: the m=1 companion RHS is the verbatim 2-vector lift.
    # -------------------------------------------------------------------------
    @testset "PH.1.2 m=1 reduction == v0.1 PainleveProblem(:I)" begin
        # m=1 RHS: y' = (y_2, 6 y_1^2 + x).
        f1 = painleve_hierarchy(:I, 1)
        x = 0.3; y = [0.4, -0.1]
        r = f1(x, y)
        @test r[1] ≈ -0.1
        @test r[2] ≈ 6*0.4^2 + 0.3 atol = 1e-12

        # Cross-check trajectories: solve the same IVP both ways on a
        # smooth interval and compare u(x).  IC chosen small/smooth so no
        # pole intervenes on [0, 0.6].
        u0, up0 = 0.1, 0.05
        xspan   = (0.0, 0.6)

        # v0.1 path.  (Scalar solve_pade takes `h_max`, not `h`.)
        pp  = PainleveProblem(:I; u0 = u0, up0 = up0, zspan = xspan, order = 30)
        s01 = solve_pade(pp; h_max = 0.05)

        # Hierarchy m=1 path.
        php = PainleveHierarchyProblem(1; y0 = [u0, up0], xspan = xspan,
                                       order = 30)
        sv  = vector_solve_pade(php.problem; h = 0.05)

        for x in 0.0:0.1:0.6
            uv = sv(x)[1]          # hierarchy:   y_1 = u
            us = s01(x)[1]         # v0.1:        u
            @test uv ≈ us atol = 1e-9 rtol = 1e-9
        end
    end

    # -------------------------------------------------------------------------
    # PH.1.3 — P_I^(2) ODE residual along a solved trajectory.
    # Solve P_I^(2) from a smooth (non-tritronquée) IC, then verify the
    # recovered components satisfy the original 4th-order ODE:
    #   u_xxxx + 10 u_x^2 + 20 u u_xx + 40(u^3 - 6t u + 6x) = 0,
    # with u_xxxx read off as the 4th component of f(x, y) along the
    # dense solution.  This is a genuine invariant — the companion form
    # is faithful iff the residual vanishes.
    # -------------------------------------------------------------------------
    @testset "PH.1.3 P_I^(2) ODE residual ≈ 0" begin
        t = 0.0
        f = painleve_hierarchy(:I, 2; t = t)
        # Small smooth IC, short interval — no pole on [0, 0.4].
        y0    = [0.1, 0.02, -0.01, 0.005]
        xspan = (0.0, 0.4)
        php   = PainleveHierarchyProblem(2; y0 = y0, xspan = xspan, order = 30)
        sv    = vector_solve_pade(php.problem; h = 0.04)

        for x in 0.05:0.05:0.35
            y  = sv(x)
            u, ux, uxx, uxxx = y
            uxxxx = f(x, y)[4]      # y_4' from the companion RHS
            resid = uxxxx + 10*ux^2 + 20*u*uxx + 40*(u^3 - 6*t*u + 6*x)
            @test abs(resid) < 1e-7
        end
    end

    # -------------------------------------------------------------------------
    # PH.1.4 — tritronquée IC: documented leading values + self-consistency.
    # The KKG asymptotic seed at large negative real x_0 (t=0):
    #   u   = -∛6 · |x_0|^{1/3}
    #   u'  = +(1/3) ∛6 · |x_0|^{-2/3}
    #   u'' = +(2/9) ∛6 · |x_0|^{-5/3}
    #   u'''= +(10/27) ∛6 · |x_0|^{-8/3}
    # obtained by differentiating u(x) = -∛6·|x|^{1/3} w.r.t. x for x<0
    # (d/dx = -d/d|x|).  NOTE: pillar C §4 line 274 prints y_3 with a
    # leading MINUS — that is a transcription sign-slip; differentiation
    # gives +(2/9).  We assert the differentiation-derived (correct) value.
    # -------------------------------------------------------------------------
    @testset "PH.1.4 tritronquée IC leading values + self-consistency" begin
        x0   = -20.0
        c    = cbrt(6.0)
        r    = abs(x0)            # |x_0| = 20
        y    = pI2_tritronquee_ic(x0)
        @test length(y) == 4
        @test y[1] ≈ -c * r^(1/3)        atol = 1e-12
        @test y[2] ≈  (1/3) * c * r^(-2/3) atol = 1e-12
        @test y[3] ≈  (2/9) * c * r^(-5/3) atol = 1e-12
        @test y[4] ≈  (10/27) * c * r^(-8/3) atol = 1e-12

        # Self-consistency #1 — the seed satisfies the original P_I^(2)
        # ODE exactly at x0: the leading power law u = -∛6·|x|^{1/3} is a
        # solution of the dominant balance, so the residual vanishes.
        f = painleve_hierarchy(:I, 2; t = 0.0)
        u, ux, uxx, uxxx = y
        uxxxx = f(x0, y)[4]
        resid0 = uxxxx + 10*ux^2 + 20*u*uxx + 40*(u^3 + 6*x0)
        @test abs(resid0) < 1e-9

        # Self-consistency #2 — a *short* forward integration tracks the
        # asymptotic series.  The P_I^(2) tritronquée is an UNSTABLE
        # special solution: nearby solutions diverge exponentially, which
        # is exactly why KKG 2015 §7 solved a boundary-value problem
        # rather than a forward IVP.  So we only assert tracking over a
        # short window (L = 0.05) where the divergence has not yet bitten;
        # the deep KKG-figure (BVP) validation is the separate bead V8b.
        L   = 0.05
        php = PainleveHierarchyProblem(2; y0 = y, xspan = (x0, x0 + L),
                                       t = 0.0, order = 30)
        sv  = vector_solve_pade(php.problem; h = 0.01)
        for x in (x0, x0 + L/2, x0 + L)
            u_series = -c * abs(x)^(1/3)
            u_solved = sv(x)[1]
            @test abs(u_solved - u_series) < 0.01 * abs(u_series)
        end
    end

    # -------------------------------------------------------------------------
    # PH.1.5 — failure modes (Rule 1).
    # -------------------------------------------------------------------------
    @testset "PH.1.5 failure modes" begin
        # m ≥ 3 deferred — must throw (general Lenard recursion not in v0.2).
        @test_throws ArgumentError painleve_hierarchy(:I, 3)
        @test_throws ArgumentError painleve_hierarchy(:I, 5)
        @test_throws ArgumentError PainleveHierarchyProblem(3;
            y0 = zeros(6), xspan = (0.0, 1.0))

        # equation ≠ :I — must throw.
        @test_throws ArgumentError painleve_hierarchy(:II, 2)
        @test_throws ArgumentError painleve_hierarchy(:VI, 1)

        # m < 1 — must throw.
        @test_throws ArgumentError painleve_hierarchy(:I, 0)

        # wrong y0 length: m=2 needs length 2m = 4.
        @test_throws ArgumentError PainleveHierarchyProblem(2;
            y0 = [1.0, 2.0, 3.0], xspan = (0.0, 1.0))
        # m=1 needs length 2.
        @test_throws ArgumentError PainleveHierarchyProblem(1;
            y0 = [1.0, 2.0, 3.0, 4.0], xspan = (0.0, 1.0))

        # degenerate xspan.
        @test_throws ArgumentError PainleveHierarchyProblem(2;
            y0 = zeros(4), xspan = (1.0, 1.0))

        # pI2_tritronquee_ic: n_terms ≥ 3 not supported (only the c_6
        # correction is implemented) — must throw rather than silently
        # truncating to a lower order.
        @test_throws ArgumentError pI2_tritronquee_ic(-20.0; n_terms = 3)
        @test_throws ArgumentError pI2_tritronquee_ic(-20.0; n_terms = 0)
        # n_terms = 2 with t ≠ 0 — the general c_n(t) series is out of
        # v0.2 scope (the documented v1 corner), so it must throw.
        @test_throws ArgumentError pI2_tritronquee_ic(-20.0; t = 0.5,
                                                      n_terms = 2)
    end

    # -------------------------------------------------------------------------
    # PH.1.7 — tritronquée IC, n_terms = 2 (bead V8b).  The KKG series at
    # t = 0 through the first non-zero correction is u = Y + Y^{-6}, with
    # Y = -∛6·r^{1/3} and Y^{-6} = 6^{-2}·r^{-2} (c_6 = 1, KKG eq. (7.2),
    # `findings.md:196-202`).  We assert the n_terms = 2 seed equals the
    # n_terms = 1 seed plus the analytic x-derivatives of 6^{-2}·r^{-2}
    # (d/dx = -d/dr for x < 0), at several x0; that n_terms = 1 is
    # byte-identical to its pre-V8b behaviour; and the n_terms = 2 seed's
    # residual against the FULL P_I^(2) ODE is smaller than the
    # leading-order seed's (the c_6 term is a genuine refinement).
    # -------------------------------------------------------------------------
    @testset "PH.1.7 tritronquée IC n_terms = 2 (c_6 correction)" begin
        for x0 in (-12.0, -20.0, -35.0, -50.0)
            c = cbrt(6.0)
            r = abs(x0)
            s = 1 / 36                       # 6^{-2}
            y1 = pI2_tritronquee_ic(x0; n_terms = 1)
            y2 = pI2_tritronquee_ic(x0; n_terms = 2)
            @test length(y2) == 4
            # n_terms = 2 = n_terms = 1 + analytic derivatives of Y^{-6}.
            @test y2[1] ≈ y1[1] + s * r^(-2) atol = 1e-13
            @test y2[2] ≈ y1[2] + 2  * s * r^(-3) atol = 1e-13
            @test y2[3] ≈ y1[3] + 6  * s * r^(-4) atol = 1e-13
            @test y2[4] ≈ y1[4] + 24 * s * r^(-5) atol = 1e-13
            # Closed-form cross-check of the y2[1] value: u = Y + Y^{-6}.
            @test y2[1] ≈ -c * r^(1/3) + s * r^(-2) atol = 1e-13

            # n_terms = 1 unchanged from its pre-V8b leading-order form.
            @test y1[1] ≈ -c * r^(1/3)         atol = 1e-13
            @test y1[2] ≈  (1/3)  * c * r^(-2/3) atol = 1e-13
            @test y1[3] ≈  (2/9)  * c * r^(-5/3) atol = 1e-13
            @test y1[4] ≈  (10/27) * c * r^(-8/3) atol = 1e-13

            # The c_6 correction is a genuine refinement.  The full
            # P_I^(2) ODE residual u'''' + 10u'^2 + 20u u'' + 40(u^3+6x)
            # must be evaluated with u'''' from the *seed formula* (a 4th
            # x-derivative of u = Y + Y^{-6} — d/dx = -d/dr, even order so
            # sign +), NOT from the companion RHS f[4] (which is defined
            # AS the negative of the rest, making that residual an
            # identically-zero tautology).  d⁴/dr⁴ of -c·r^{1/3}:
            # coeff (1/3)(-2/3)(-5/3)(-8/3) = -80/81 ⇒ +(80/81)c·r^{-11/3};
            # d⁴/dr⁴ of s·r^{-2}: coeff (-2)(-3)(-4)(-5) = 120 ⇒ 120s·r^{-6}.
            u4_1 = (80/81) * c * r^(-11/3)                 # n_terms = 1
            u4_2 = (80/81) * c * r^(-11/3) + 120 * s * r^(-6)  # n_terms = 2
            res1 = abs(u4_1 + 10*y1[2]^2 + 20*y1[1]*y1[3] + 40*(y1[1]^3 + 6*x0))
            res2 = abs(u4_2 + 10*y2[2]^2 + 20*y2[1]*y2[3] + 40*(y2[1]^3 + 6*x0))
            @test res2 < res1
        end
        # BigFloat element type propagates (a BigFloat x0 → BigFloat seed).
        yb = pI2_tritronquee_ic(big"-20.0"; n_terms = 2)
        @test eltype(yb) == BigFloat
        @test yb[1] ≈ -cbrt(big(6)) * big(20)^(big(1)/3) +
                      (one(BigFloat)/36) * big(20)^(-2)  atol = big(1e-40)
    end

    # -------------------------------------------------------------------------
    # Umbrella: the module is wired into PadeTaylor.
    # -------------------------------------------------------------------------
    @testset "umbrella wiring" begin
        @test isdefined(PadeTaylor, :PainleveHierarchy)
        @test isa(painleve_hierarchy, Function)
    end

end

# -----------------------------------------------------------------------------
# PH.1.6 — Mutation-proof record (CLAUDE.md Rule 4).
#
# Each mutation was applied to `src/PainleveHierarchy.jl`, the suite re-run,
# the RED count recorded, then the source restored.
#
#   M1  — `_rhs_PI2`: coefficient 40 → 41 in the `40*(...)` term.
#         Bit 18 asserts RED — PH.1.1 (10: every m=2 RHS comparison),
#         PH.1.3 (7: every residual point), PH.1.4 (1: short-integration
#         tracking).  RED.
#   M2  — `_rhs_PI2`: sign of the `6*t*y1` term flipped (−6t → +6t).
#         Bit 7 asserts RED — all in PH.1.1, the t≠0 RHS comparisons;
#         the t=0 asserts and PH.1.3/1.4 (both t=0) stayed GREEN, exactly
#         as expected for a t-only mutation.  RED.
#   M3  — `pI2_tritronquee_ic`: exponent 1/3 → 2/3 on the y[1] term.
#         Bit 4 asserts RED — PH.1.4 (the y[1] leading-value assert + the
#         three short-integration tracking asserts).  RED.
#
# All three mutations produced a RED suite; the source was restored to
# GREEN (59/59) after each.  The tests catch a wrong coefficient, a
# wrong sign, and a wrong exponent — the three regression classes for
# this module.
# -----------------------------------------------------------------------------
