# test/fw_fig_41_test.jl — FW 2011 Fig 4.1 step-(i) BVP quantitative pin.
#
# Bead `padetaylor-0c3`.  See `docs/worklog/016-fw-fig-41-pin.md`.
#
# Reproduces the tritronquée FW eq. 4.1 reference values via the FW Fig 4.1
# step-(i) recipe: Chebyshev-Newton BVP for `u'' = 6u² + z` on the
# vertical imaginary-axis segment `[-20i, +20i]`, with the leading-term
# asymptote `u(z) = -√(-z/6)` (FW eq. 1.2, negative branch) used as
# Dirichlet BCs at both endpoints.
#
# References:
#   - `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`
#     :48 (eq. 1.1: u'' = 6u² + z), :54-55 (eq. 1.2: u ~ ±√(-z/6)),
#     :214-222 (Fig 4.1 caption: 3-step composition), :226-227
#     (FW eq. 4.1: u(0), u'(0) tritronquée reference values).
#   - `docs/figure_catalogue.md` §1 row 4.1.
#
# Why this lives in a standalone test file (not lattice_dispatcher_test.jl):
# the FW Fig 4.1 algorithm is NOT the horizontal-row BVP-fill of
# `lattice_dispatch_solve` (line 190).  It's a VERTICAL-axis BVP whose
# output FEEDS pole-field IVP integrators (step (ii)) and an interior
# BVP fill (step (iii)).  v1 ships only the step-(i) pin; the composition
# of steps (ii) and (iii) is downstream work (no bead filed; the lattice
# dispatcher v1 covers the line-190 pattern, and step (iii) is structurally
# similar).

using Test
using PadeTaylor

@testset "FW 2011 Fig 4.1 step-(i) BVP: tritronquée pin via leading-term BCs" begin

    # Painlevé I in BVP form: u'' = F(z, u) where F(z, u) = 6u² + z.
    # ∂F/∂u = 12u (diagonal Jacobian — the BVP solver's Newton step
    # uses this directly, no automatic differentiation).
    bvp_f(z, u)      = 6 * u^2 + z
    bvp_∂f_∂u(z, u)  = 12 * u

    # Leading-term asymptote.  Negative branch is the tritronquée
    # selector (4 pole-free sectors approached from `u → 0⁻` along the
    # negative real axis ⇒ u(0) ≈ -0.188 < 0, consistent with -√).
    # FW2011...md:54.  +√ branch diverges in Newton (verified empirically;
    # see worklog 016 §"BC sign drives basin selection").
    leading(z) = -sqrt(-z / 6)

    z_a = -20.0im
    z_b = +20.0im
    u_a = leading(z_a)
    u_b = leading(z_b)

    # N=240 calibrated empirically (worklog 016 §"Spectral convergence
    # probe"): err_u(0) ~ 3.5e-13, err_up(0) ~ 5.3e-11.  Both ≤ 1e-10.
    # Wall ≈ 100 ms for the single BVP solve.  Smaller N undershoots
    # the pin (N=220: err_up ≈ 4.9e-10 — fails).
    sol = bvp_solve(bvp_f, bvp_∂f_∂u, z_a, z_b, u_a, u_b;
                     N = 240, initial_guess = leading,
                     tol = 1e-13, maxiter = 20)

    @testset "FF.1.1: Newton converged" begin
        @test sol.iterations ≤ 10
        @test isfinite(sol.residual_inf)
        # `residual_inf` reflects the truncation-error floor at the
        # converged Newton point, NOT the Newton step norm.  At N=240
        # it sits around 1e-8 — well above the `eps^(3/4)` step
        # tolerance used by the BVP's step-norm Newton criterion
        # (worklog 006 §"Step-norm Newton").  We don't gate on residual.
    end

    @testset "FF.1.2: u(0), u'(0) match FW eq. 4.1 to ≤ 1e-10" begin
        # FW eq. 4.1 — tritronquée ICs at z=0, to 16 digits (FW2011...md:226).
        u_at_0_FW  = -0.1875543083404949
        up_at_0_FW =  0.3049055602612289

        u_0, up_0 = sol(0.0 + 0.0im)

        @test abs(real(u_0)  - u_at_0_FW)  ≤ 1e-10
        @test abs(real(up_0) - up_at_0_FW) ≤ 1e-10

        # Schwarz reflection (real-coefficient ODE + conjugate-symmetric
        # BCs ⇒ u, u' real on the real axis): imag parts negligible.
        @test abs(imag(u_0))  ≤ 1e-10
        @test abs(imag(up_0)) ≤ 1e-10
    end

    @testset "FF.1.3: u(+20i) = imposed Dirichlet BC (structural consistency)" begin
        # Sanity: the barycentric eval at the segment endpoint must
        # return the imposed BC to within roundoff.  The exact Dirichlet
        # enforcement is a property of the Chebyshev collocation, not a
        # numerical assertion against an independent oracle.
        u_top, _ = sol(z_b)
        @test abs(u_top - u_b) ≤ 1e-13
    end

    @testset "FF.1.4: u'(+20i) sanity (finite + leading-term magnitude)" begin
        # Pinning u'(20i) to an independent oracle is deferred — the BVP
        # imposes a Dirichlet BC on u alone, so u'(20i) is whatever the
        # converged solution returns.  Per FW2011...md:229, the leading-
        # term BC is "crude" and the BVP "smooths out" the BC error over
        # the long segment to give u(0)/u'(0) accurate to ≤1e-20; that
        # smoothing implies u'(20i) does NOT exactly equal the leading-
        # term derivative.  We assert structural sanity only:
        #   - leading-term derivative at z=20i:  d/dz[-√(-z/6)] = 1/(12·√(-z/6))
        #     ≈ 0.0323 + 0.0323i (verified by hand; worklog 016).
        #   - The BVP gives u'(20i) ≈ 0.0324 + 0.0325i — same order of
        #     magnitude, ~2.5e-4 absolute difference (the o(1) correction).
        _, up_top = sol(z_b)
        up_leading = 1 / (12 * sqrt(-z_b / 6))
        @test isfinite(real(up_top)) && isfinite(imag(up_top))
        @test abs(up_top - up_leading) ≤ 1e-3   # o(1) correction; not 1e-10
        @test 0.01 < abs(up_top) < 0.1
    end

    @testset "FF.1.5: Schwarz symmetry u(-iy) = conj(u(iy)) on imag axis" begin
        # The ODE preserves complex conjugation (real coefficients) and
        # the BCs are conjugate-symmetric (u_a = conj(u_b)).  The converged
        # BVP solution must therefore satisfy u(-z) = conj(u(z)) on the
        # imaginary-axis segment.  Tolerance is "spectral-truncation-floor"
        # × small, dominated by barycentric eval roundoff.
        max_asym = 0.0
        for y in (2.0, 5.0, 10.0, 15.0, 19.0)
            u_up, _ = sol(+im * y)
            u_dn, _ = sol(-im * y)
            max_asym = max(max_asym, abs(u_up - conj(u_dn)))
        end
        @test max_asym ≤ 1e-12
    end

    @testset "FF.1.6: +√ branch — wrong basin, Newton diverges" begin
        # Mutation-proof companion (verified manually 2026-05-13): the
        # positive-square-root BC `leading_pos(z) = +√(-z/6)` does NOT
        # select the tritronquée basin.  At z = 0 the leading term
        # vanishes from BELOW (-0.188 < 0), so the negative branch is
        # the one consistent with the tritronquée's near-origin shape.
        # The +√ BC drives Newton into a divergent basin (‖Δu‖ grows;
        # ‖R‖ → 1e4+).  This test pins that fact.
        leading_pos(z) = +sqrt(-z / 6)
        @test_throws ErrorException bvp_solve(
            bvp_f, bvp_∂f_∂u, z_a, z_b, leading_pos(z_a), leading_pos(z_b);
            N = 240, initial_guess = leading_pos,
            tol = 1e-13, maxiter = 10)
    end

end # @testset FW 2011 Fig 4.1
