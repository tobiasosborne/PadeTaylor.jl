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
    PH.1.4  tritronquée IC     — pI2_tritronquee_ic returns the CORRECT-
                                  SIGN leading values (u > 0 on x < 0);
                                  the seed's P_I^(2) ODE residual is small
                                  and below the leading-order bound 6/|x|
    PH.1.5  failure modes      — m ≥ 3 throws; equation ≠ :I throws;
                                  wrong y0 length throws; n_terms ∉ {1,2}
                                  throws; n_terms = 2 with t ≠ 0 throws;
                                  abs(x0) < 1 throws
    PH.1.7  tritronquée IC n=2 — pI2_tritronquee_ic(·; n_terms = 2) at
                                  t = 0 equals Y + Y^{-6} (c_6 = 1) and
                                  its analytic derivatives; the n=2 seed's
                                  ODE residual is below 20/|x|^4 — orders
                                  of magnitude tighter than n=1 (bead V8b)
    PH.1.8  tritronquée IC      — pI2_tritronquee_ic generalised to any
            off-axis             complex x0 on the V_0 sheet (bead 0ln.32);
                                  off-axis rays φ ∈ {2π/3,π,4π/3} give a
                                  small, R-decaying P_I^(2) ODE residual,
                                  and the φ=π ray reproduces the negative-
                                  real-axis seed
    PH.1.6  mutation-proof     — recorded at end of file

The decisive ground truth for the tritronquée seed (PH.1.4, PH.1.7) is
the P_I^(2) ODE itself.  The leading dominant balance 40(u³ + 6x) = 0
gives u³ = -6x; on the negative real axis (x < 0) u³ = 6|x| > 0, so the
REAL root is u = +∛6·|x|^{1/3} — POSITIVE.  The seed (u, u', u'', u''')
is plugged into u'''' + 10u'^2 + 20u·u'' + 40(u³ + 6x) with u'''' taken
from the seed's CLOSED FORM (a genuine residual, not the f[4] tautology)
and the residual asserted against its leading-order decay rate.

Oracles: the hand-written companion RHS (PH.1.1); v0.1
`PainleveProblem(:I)` trajectory (PH.1.2); the original ODE residual
computed from the solved jet (PH.1.3); the P_I^(2) ODE residual of the
asymptotic seed, with u'''' from the seed's closed form (PH.1.4, PH.1.7).

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
    # recovered components satisfy the original 4th-order ODE
    #   u_xxxx + 10 u_x^2 + 20 u u_xx + 40(u^3 - 6t u + 6x) = 0.
    #
    # bead-0ln.30 caveat: u_xxxx must NOT be read off as f(x,y)[4] of the
    # solved state.  f[4] is DEFINED as -(10y2^2 + 20y1y3 + 40(...)), so
    # the residual f[4] + 10y2^2 + 20y1y3 + 40(...) ≡ 0 for ANY y — a
    # tautology that says nothing about whether `sv` actually solved the
    # ODE.  The genuine invariant tested here instead is the companion-
    # CHAIN consistency: along the solved trajectory the dense components
    # must satisfy y_{k+1} = d/dx y_k.  We finite-difference the solved
    # third component y_4 = u''' to get an INDEPENDENT u'''' and plug
    # THAT into the ODE — a residual that goes large the moment the
    # integrator returns an inconsistent jet.
    # -------------------------------------------------------------------------
    @testset "PH.1.3 P_I^(2) ODE residual ≈ 0" begin
        t = 0.0
        # Small smooth IC, short interval — no pole on [0, 0.4].
        y0    = [0.1, 0.02, -0.01, 0.005]
        xspan = (0.0, 0.4)
        php   = PainleveHierarchyProblem(2; y0 = y0, xspan = xspan, order = 30)
        sv    = vector_solve_pade(php.problem; h = 0.04)

        for x in 0.05:0.05:0.35
            y  = sv(x)
            u, ux, uxx, uxxx = y
            # INDEPENDENT u'''' — central-difference the solved u''' (the
            # 4th companion component); NOT the companion RHS f[4].
            dx     = 1e-4
            uxxxx  = (sv(x + dx)[4] - sv(x - dx)[4]) / (2*dx)
            resid  = uxxxx + 10*ux^2 + 20*u*uxx + 40*(u^3 - 6*t*u + 6*x)
            @test abs(resid) < 1e-5
            # Companion-chain consistency: y_{k+1} = d/dx y_k along sv.
            @test (sv(x + dx)[1] - sv(x - dx)[1]) / (2*dx) ≈ ux  atol = 1e-6
            @test (sv(x + dx)[2] - sv(x - dx)[2]) / (2*dx) ≈ uxx atol = 1e-6
        end
    end

    # -------------------------------------------------------------------------
    # PH.1.4 — tritronquée IC: CORRECT-SIGN leading values + genuine ODE
    # residual.  The KKG asymptotic seed at large negative real x_0 (t=0)
    # is u = +∛6·|x_0|^{1/3} — POSITIVE on x < 0.  Ground truth: the
    # P_I^(2) leading balance 40(u³ + 6x) = 0 gives u³ = -6x = 6|x| > 0,
    # so the real cube root is positive (bead 0ln.31; KKG's uniform
    # u ~ -∛6·x^{1/3} with the *real* cube root, x^{1/3} = -|x|^{1/3} for
    # x<0).  Differentiating u = +∛6·|x|^{1/3} for x<0 (d/dx = -d/d|x|)
    # flips every odd-order derivative:
    #   u   = +∛6 · |x_0|^{1/3}
    #   u'  = -(1/3) ∛6 · |x_0|^{-2/3}
    #   u'' = -(2/9) ∛6 · |x_0|^{-5/3}
    #   u'''= -(10/27) ∛6 · |x_0|^{-8/3}
    #
    # The load-bearing assertion is a GROUND-TRUTH ODE-residual check
    # (Rule 5).  Plug (u,u',u'',u''') into the original 4th-order ODE
    #   u'''' + 10u'^2 + 20u·u'' + 40(u³ - 6t·u + 6x)
    # with u'''' taken from the seed's CLOSED FORM, NOT from the companion
    # RHS f(x,y)[4].  This is the bead-0ln.30 fix: f[4] is DEFINED as the
    # negative of the rest of the row, so resid = f[4] + rest ≡ 0 is an
    # identity that holds for ANY y — a tautology, not a test.  The
    # closed-form u'''' makes it a genuine, non-trivial residual.
    # d⁴/dx⁴ of u = +∛6·|x|^{1/3} (x<0, d/dx=-d/d|x|, 4th order ⇒ even ⇒
    # sign +): coeff (1/3)(-2/3)(-5/3)(-8/3) = -80/81 ⇒ u'''' =
    # -(80/81)·∛6·|x|^{-11/3}.  The leading-order seed is exact only at
    # dominant-balance order, so the residual is the leading truncation
    # error; empirically |resid| ≈ 5/|x| (decaying), so 6/|x| is a safe
    # pin (verified by /tmp/verify_scale.jl: 4.06 at x=-20, 2.99 at -50).
    # -------------------------------------------------------------------------
    @testset "PH.1.4 tritronquée IC: correct sign + genuine ODE residual" begin
        c = cbrt(6.0)
        for x0 in (-12.0, -20.0, -35.0, -50.0)
            r = abs(x0)
            y = pI2_tritronquee_ic(x0)
            @test length(y) == 4
            # Correct-sign leading values: u>0, u'<0, u''<0, u'''<0.
            @test y[1] ≈  c * r^(1/3)          atol = 1e-12
            @test y[2] ≈ -(1/3)  * c * r^(-2/3) atol = 1e-12
            @test y[3] ≈ -(2/9)  * c * r^(-5/3) atol = 1e-12
            @test y[4] ≈ -(10/27) * c * r^(-8/3) atol = 1e-12
            @test y[1] > 0          # u is POSITIVE on the negative axis

            # Genuine ODE residual — u'''' from the seed's closed form.
            u, ux, uxx, uxxx = y
            uxxxx = -(80/81) * c * r^(-11/3)
            resid = uxxxx + 10*ux^2 + 20*u*uxx + 40*(u^3 + 6*x0)
            @test abs(resid) < 6 / r        # leading-order truncation bound
            # The OLD (negative) branch is NOT a solution: its residual is
            # ≈ -480·|x| ~ thousands.  Pin that the correct seed beats it
            # by a wide margin (a few orders of magnitude).
            @test abs(resid) < 1.0
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
        # pI2_tritronquee_ic: |x0| < 1 must throw — the KKG asymptotic
        # series u ~ -∛6·x^{1/3} + … is meaningless at small |x0| (the
        # V8/V8b figures seed at |x0| ≈ 20–50); fail-fast (Rule 1).
        @test_throws ArgumentError pI2_tritronquee_ic(0.0)
        @test_throws ArgumentError pI2_tritronquee_ic(-0.5)
        @test_throws ArgumentError pI2_tritronquee_ic(0.3im)
        @test_throws ArgumentError pI2_tritronquee_ic(-0.5; n_terms = 2)
    end

    # -------------------------------------------------------------------------
    # PH.1.7 — tritronquée IC, n_terms = 2 (bead V8b).  The KKG series at
    # t = 0 through the first non-zero correction is u = Y + Y^{-6}, with
    # Y = +∛6·r^{1/3} (the CORRECT-SIGN real branch on x<0, bead 0ln.31)
    # and Y^{-6} = 6^{-2}·r^{-2} (c_6 = 1, KKG eq. (7.2),
    # `findings.md:196-202`).  The c_6 correction Y^{-6} is an EVEN power
    # of Y, hence branch-insensitive — the sign fix touches only the
    # n_terms = 1 leading block; the Δ-pieces are unchanged.
    #
    # We assert the n_terms = 2 seed equals the n_terms = 1 seed plus the
    # analytic x-derivatives of 6^{-2}·r^{-2} (d/dx = -d/dr for x < 0);
    # the correct-sign leading values; and — the load-bearing GROUND-
    # TRUTH check (Rule 5) — that the n_terms = 2 seed's P_I^(2) ODE
    # residual is below 20/|x|^4, orders of magnitude tighter than the
    # n_terms = 1 seed's ≈ 5/|x|.  The c_6 term IS a genuine refinement.
    # u'''' is taken from the seed's CLOSED FORM (bead 0ln.30 — never the
    # companion f[4], which gives an identically-zero tautology):
    #   d⁴/dx⁴ of +c·r^{1/3}: coeff (1/3)(-2/3)(-5/3)(-8/3) = -80/81 ⇒
    #     -(80/81)c·r^{-11/3} (even order, d/dx=-d/dr sign cancels);
    #   d⁴/dx⁴ of s·r^{-2}: coeff (-2)(-3)(-4)(-5) = 120 ⇒ 120s·r^{-6}.
    # Empirical n=2 residual (verify_scale.jl): |R2|·|x|^4 ∈ [8,13], so
    # 20/|x|^4 is a safe pin; n=1: |R1|·|x| ∈ [3,5], pin 6/|x|.
    # -------------------------------------------------------------------------
    @testset "PH.1.7 tritronquée IC n_terms = 2 (c_6 correction)" begin
        for x0 in (-12.0, -20.0, -35.0, -50.0)
            c = cbrt(6.0)
            r = abs(x0)
            s = 1 / 36                       # 6^{-2}
            y1 = pI2_tritronquee_ic(x0; n_terms = 1)
            y2 = pI2_tritronquee_ic(x0; n_terms = 2)
            @test length(y2) == 4
            # n_terms = 2 = n_terms = 1 + analytic derivatives of Y^{-6}
            # (branch-insensitive: Δ-pieces all positive, unchanged).
            @test y2[1] ≈ y1[1] + s * r^(-2) atol = 1e-13
            @test y2[2] ≈ y1[2] + 2  * s * r^(-3) atol = 1e-13
            @test y2[3] ≈ y1[3] + 6  * s * r^(-4) atol = 1e-13
            @test y2[4] ≈ y1[4] + 24 * s * r^(-5) atol = 1e-13
            # Closed-form cross-check of y2[1]: u = Y + Y^{-6}, Y>0.
            @test y2[1] ≈  c * r^(1/3) + s * r^(-2) atol = 1e-13

            # n_terms = 1 correct-sign leading-order form: u>0, u'<0,…
            @test y1[1] ≈  c * r^(1/3)          atol = 1e-13
            @test y1[2] ≈ -(1/3)  * c * r^(-2/3) atol = 1e-13
            @test y1[3] ≈ -(2/9)  * c * r^(-5/3) atol = 1e-13
            @test y1[4] ≈ -(10/27) * c * r^(-8/3) atol = 1e-13
            @test y1[1] > 0 && y2[1] > 0

            # GROUND-TRUTH ODE residual — u'''' from the seed's closed
            # form (NOT companion f[4]: that residual is identically 0).
            u4_1 = -(80/81) * c * r^(-11/3)                 # n_terms = 1
            u4_2 = -(80/81) * c * r^(-11/3) + 120 * s * r^(-6)  # n_terms = 2
            res1 = abs(u4_1 + 10*y1[2]^2 + 20*y1[1]*y1[3] + 40*(y1[1]^3 + 6*x0))
            res2 = abs(u4_2 + 10*y2[2]^2 + 20*y2[1]*y2[3] + 40*(y2[1]^3 + 6*x0))
            # The c_6 term is a genuine refinement: n=2 residual is far
            # smaller than n=1, and below the documented |x|^{-4} bound.
            @test res1 < 6  / r           # n=1 leading-order bound
            @test res2 < 20 / r^4         # n=2 c_6-corrected bound
            @test res2 < res1
        end
        # BigFloat element type propagates (a BigFloat x0 → BigFloat seed).
        yb = pI2_tritronquee_ic(big"-20.0"; n_terms = 2)
        @test eltype(yb) == BigFloat
        @test yb[1] ≈ cbrt(big(6)) * big(20)^(big(1)/3) +
                      (one(BigFloat)/36) * big(20)^(-2)  atol = big(1e-40)
    end

    # -------------------------------------------------------------------------
    # PH.1.8 — tritronquée IC on OFF-AXIS rays of the V_0 sheet (bead
    # 0ln.32).  pI2_tritronquee_ic is generalised from negative-real-axis-
    # only to any complex x0 on the V_0 sheet.  The branch places the
    # negative real axis at sheet angle θ = 3π: from the principal angle
    # φ₀ = angle(x0) ∈ (-π,π], θ = φ₀≤0 ? φ₀+4π : φ₀+2π ∈ (2π,4π], and
    # Y = -∛6·|x0|^{1/3}·exp(iθ/3).  The higher components are the branch-
    # free rational functions y₂ = Y/(3x0), y₃ = -2Y/(9x0²),
    # y₄ = 10Y/(27x0³); the c_6 correction adds (x0^{-2}, -2x0^{-3},
    # 6x0^{-4}, -24x0^{-5})/36 (an even power of Y, single-valued in x0).
    #
    # Two ground-truth checks (Rule 5):
    #   (a) the three test rays φ ∈ {2π/3, π, 4π/3} map to θ ∈ {8π/3, 3π,
    #       10π/3} — all inside the V_0 pole-free sector
    #       θ ∈ (3π − 6π/7+…, 3π + 6π/7−…) — and the φ = π ray reproduces
    #       the negative-real-axis seed (continuity at the axis);
    #   (b) the load-bearing assertion: the seed plugged into the P_I^(2)
    #       ODE u'''' + 10u'² + 20u·u'' + 40(u³ + 6x) has a SMALL residual
    #       that decays with R, with u'''' from the seed's CLOSED FORM
    #       (NOT companion f[4] — that is the 0ln.30 tautology).
    #       u'''' of Y = -∛6·x^{1/3}: coeff (1/3)(-2/3)(-5/3)(-8/3) = -80/81
    #       and Y = -∛6·x^{1/3} ⇒ u4_Y = (80/81)·∛6·x0^{1/3}·x0^{-4};
    #       u'''' of x0^{-2}/36: coeff (-2)(-3)(-4)(-5) = 120 ⇒ (120/36)x0^{-6}.
    #   Empirically (verified in-session): |R1|·R ∈ [3.5,4.5] for n=1 and
    #   |R2|·R⁴ ∈ [8.7,11] for n=2 — the SAME decay as on the negative
    #   real axis (PH.1.4/1.7), so the bounds 6/R and 20/R⁴ carry over.
    # -------------------------------------------------------------------------
    @testset "PH.1.8 tritronquée IC on off-axis V_0-sheet rays" begin
        c = cbrt(6.0)
        for R in (15.0, 30.0), φ in (2π/3, 1.0π, 4π/3)
            x0 = R * exp(im * φ)
            for nt in (1, 2)
                y = pI2_tritronquee_ic(x0; n_terms = nt)
                @test length(y) == 4
                @test eltype(y) <: Complex

                # Reconstruct the branch independently and cross-check the
                # leading component Y and the rational derivatives, plus
                # the branch-free c_6 correction block when nt == 2.
                φ0 = angle(x0)
                θ  = φ0 ≤ 0 ? φ0 + 4π : φ0 + 2π
                Y  = -c * abs(x0)^(1/3) * exp(im * θ / 3)
                δ  = nt == 2 ?
                     [x0^(-2), -2x0^(-3), 6x0^(-4), -24x0^(-5)] ./ 36 :
                     zeros(ComplexF64, 4)
                @test y[1] ≈ Y            + δ[1] atol = 1e-12
                @test y[2] ≈ Y/(3x0)      + δ[2] atol = 1e-12
                @test y[3] ≈ -2Y/(9x0^2)  + δ[3] atol = 1e-12
                @test y[4] ≈ 10Y/(27x0^3) + δ[4] atol = 1e-12

                # GROUND-TRUTH ODE residual — u'''' from the closed form.
                u, ux, uxx, _ = y
                u4 = (80/81) * c * abs(x0)^(1/3) * exp(im*θ/3) * x0^(-4)
                nt == 2 && (u4 += (120/36) * x0^(-6))
                resid = u4 + 10*ux^2 + 20*u*uxx + 40*(u^3 + 6*x0)
                if nt == 1
                    @test abs(resid) < 6 / R          # leading-order bound
                else
                    @test abs(resid) < 20 / R^4       # c_6-corrected bound
                end
            end
            # n=2 residual must beat n=1 (the c_6 term is a refinement).
            y1 = pI2_tritronquee_ic(x0; n_terms = 1)
            y2 = pI2_tritronquee_ic(x0; n_terms = 2)
            φ0 = angle(x0); θ = φ0 ≤ 0 ? φ0 + 4π : φ0 + 2π
            r1 = abs((80/81)*c*abs(x0)^(1/3)*exp(im*θ/3)*x0^(-4) +
                     10*y1[2]^2 + 20*y1[1]*y1[3] + 40*(y1[1]^3 + 6*x0))
            r2 = abs((80/81)*c*abs(x0)^(1/3)*exp(im*θ/3)*x0^(-4) +
                     (120/36)*x0^(-6) +
                     10*y2[2]^2 + 20*y2[1]*y2[3] + 40*(y2[1]^3 + 6*x0))
            @test r2 < r1
        end

        # Continuity at the negative real axis: the φ = π ray (a complex
        # x0 = R·e^{iπ}) seed agrees with the real negative-axis seed.
        for R in (15.0, 30.0)
            y_axis = pI2_tritronquee_ic(-R; n_terms = 2)        # real path
            y_ray  = pI2_tritronquee_ic(R*exp(im*1.0π); n_terms = 2)
            @test maximum(abs.(y_ray .- y_axis)) < 1e-12
        end

        # Element type: Complex{BigFloat} x0 → Complex{BigFloat} seed.
        yb = pI2_tritronquee_ic(big(20.0) * exp(im * big(2π/3)); n_terms = 1)
        @test eltype(yb) == Complex{BigFloat}
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
#         tracking).  RED.  [pre-0ln.31 record; PH.1.3/1.4 since reshaped]
#   M2  — `_rhs_PI2`: sign of the `6*t*y1` term flipped (−6t → +6t).
#         Bit 7 asserts RED — all in PH.1.1, the t≠0 RHS comparisons;
#         the t=0 asserts and PH.1.3/1.4 (both t=0) stayed GREEN, exactly
#         as expected for a t-only mutation.  RED.
#   M3  — `pI2_tritronquee_ic`: exponent 1/3 → 2/3 on the y[1] term.
#         Bit 4 asserts RED — PH.1.4 (the y[1] leading-value assert + the
#         three short-integration tracking asserts).  RED.  [pre-0ln.31]
#   M4  — bead 0ln.31 sign-fix mutation.  `pI2_tritronquee_ic`: the four
#         corrected leading-order seed signs flipped back to the OLD
#         (wrong) branch — y[1] = -c·r^{1/3}, y[2..4] sign-flipped — the
#         exact bug 0ln.31 fixes.  RESULT: 61 of 159 assertions RED,
#         dominated by the PH.1.4 and PH.1.7 GROUND-TRUTH ODE-residual
#         checks: `abs(resid) < 6/r` evaluated to 9600.2 < 0.3 at x=-20
#         and 5760.4 < 0.5 at x=-12 — the residual jumps from ~0.2 (the
#         correct leading-order truncation error) to ~thousands (the
#         wrong branch is not a P_I^(2) solution at all, residual ≈
#         -480·|x|).  Confirms the new residual assertions are GENUINE:
#         a tautological `f[4]`-based check would have stayed GREEN under
#         this mutation (f[4] ≡ -(rest) regardless of sign).  Restored to
#         GREEN (159/159).
#   M5  — bead 0ln.32 generalisation.  `pI2_tritronquee_ic` general
#         (off-axis) branch: the V_0-sheet leading order
#         `Y = -∛6·|z|^{1/3}·exp(iθ/3)` mutated to `exp(iθ)` — the
#         cube-root exponent dropped, so Y is no longer a cube root of
#         -6x at all.  RESULT: 40 of 253 assertions RED, the PH.1.8
#         off-axis block: the GROUND-TRUTH ODE residual `abs(resid) <
#         6/R` jumped to 6235.6 (R=15) / 12470.9 (R=30) on the φ=2π/3
#         and 4π/3 rays — the same ≈ -480·|x| explosion M4 exposed on
#         the real axis.  (The φ=π ray stayed GREEN: there exp(iθ) =
#         exp(i3π) = -1 = exp(iπ) coincidentally, the negative-real-axis
#         special case.)  Confirms the off-axis residual check is
#         GENUINE.  Restored to GREEN (253/253).
#         FALSE START — a first mutation, dropping the `+2π`/`+4π`
#         sheet shift (`θ = φ₀`), did NOT bite the ODE-residual
#         assertion: `θ → θ − 2π` multiplies Y by `exp(-2πi/3)`, a cube
#         root of unity, so the mutated Y is still *a* cube root of -6x
#         and the leading dominant balance still holds — the residual
#         stayed ≈ 0.3.  It only bit the explicit branch cross-checks
#         (y[1..4] ≈ Y…).  Lesson: the dominant balance has a 3-fold
#         cube-root ambiguity the residual check cannot see; the
#         load-bearing mutation must perturb the *power* (θ/3 → θ), not
#         the sheet index.
#
# All five mutations produced a RED suite; the source was restored to
# GREEN after each.  The tests catch a wrong coefficient (M1), a wrong
# t-sign (M2), a wrong exponent (M3), the wrong tritronquée branch on
# the real axis (M4), and a wrong off-axis cube-root power (M5) — the
# last two detectable only by the closed-form-u'''' residual check
# (not the f[4] tautology).
# -----------------------------------------------------------------------------
