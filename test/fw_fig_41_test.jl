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

# ======================================================================
# Steps (ii) + (iii) — bead `padetaylor-gky`, worklog 029.
#
# Step (i) above pins the imaginary-axis BVP.  Steps (ii) and (iii)
# compose it into the full Fig 4.1 (`figures/fw2011_fig_4_1.jl`):
#   (ii)  run out the pole field with `edge_gated_pole_field_solve`
#         from the BVP-derived IC at z = 0;
#   (iii) fill the smooth band per grid line with `bvp_solve`.
#
# There is no in-tree quantitative oracle for PI tritronquée pole
# *locations* (FW Table 5.1 is for the equianharmonic ℘, not PI — see
# worklog 016 §"What is NOT shipped").  But the composition admits two
# oracle-free cross-validations and one self-consistency check, all
# pinned here:
#   FF.2.1 — the step-(ii) IVP run-out and the step-(i) BVP spine must
#            agree where they overlap (independent methods — spectral
#            BVP vs. Taylor–Padé IVP — on the same solution).
#   FF.2.2 — the step-(ii) run-out from a real IC is conjugate-
#            symmetric near z = 0 (real ODE coefficients + real IC).
#   FF.2.3 — the step-(iii) fill primitive: `bvp_solve` on a sub-
#            segment of the imaginary axis, given exact spine values
#            as BCs, reproduces the spine in the interior.
# ======================================================================
@testset "FW 2011 Fig 4.1 steps (ii)+(iii): BVP+IVP composition" begin

    pI_ivp(z, u, up) = 6 * u^2 + z
    bvp_f(z, u)      = 6 * u^2 + z
    bvp_∂f_∂u(z, u)  = 12 * u
    leading(z)       = -sqrt(-z / 6)

    # Step (i): the imaginary-axis BVP spine (the recipe pinned above).
    z_a, z_b = -20.0im, 20.0im
    spine = bvp_solve(bvp_f, bvp_∂f_∂u, z_a, z_b, leading(z_a), leading(z_b);
                      N = 240, initial_guess = leading, tol = 1e-13, maxiter = 20)
    u0, up0 = spine(0.0 + 0.0im)

    # Step (ii): edge-gated pole-field run-out from the BVP-derived IC.
    # A small lattice over [-3,3] × [-3,9] — the seed disc around z = 0
    # covers the imaginary axis near the origin (overlap with the
    # spine) and reaches the first tritronquée pole at z ≈ 2.07.
    xs = range(-3.0, 3.0; step = 0.5)
    ys = range(-3.0, 9.0; step = 0.5)
    i_axis = findfirst(x -> abs(x) < 0.25, xs)
    prob0 = PadeTaylorProblem(pI_ivp, (real(u0), real(up0)),
                              (0.0 + 0.0im, 3.0 + 0.0im); order = 30)
    egs0 = edge_gated_pole_field_solve(prob0, xs, ys; h = 0.5, grow_rings = 3)

    @testset "FF.2.1: step-(ii) IVP run-out ≡ step-(i) BVP spine" begin
        # On the imaginary axis near z = 0 both methods describe the
        # same near-tritronquée solution: the IVP propagates (u0,u'0)
        # via Taylor–Padé, the BVP solves the 2-point problem
        # spectrally.  Where the IVP run-out has a value on the x ≈ 0
        # column, it must match `spine`.
        worst = 0.0
        ncheck = 0
        for j in eachindex(ys)
            u_ivp = egs0.u_grid[i_axis, j]
            isfinite(u_ivp) || continue
            u_bvp = spine(im * ys[j])[1]
            worst = max(worst, abs(u_ivp - u_bvp))
            ncheck += 1
        end
        @test ncheck ≥ 5                       # genuine overlap exists
        @test worst ≤ 1e-6                     # independent methods agree
    end

    @testset "FF.2.2: step-(ii) run-out is conjugate-symmetric near z=0" begin
        # Real ODE coefficients + real IC ⇒ u(z̄) = conj(u(z)).  The
        # path-network's RNG-shuffled target order makes the visited
        # tree asymmetric, so this holds tightly only near the IC
        # (worklog 014); we check the seed neighbourhood |y| ≤ 2.
        worst = 0.0
        ncheck = 0
        for j in eachindex(ys), i in eachindex(xs)
            y = ys[j]
            (0 < y ≤ 2.0) || continue
            jm = findfirst(yy -> abs(yy + y) < 0.25, ys)   # row at -y
            jm === nothing && continue
            u_up = egs0.u_grid[i, j]
            u_dn = egs0.u_grid[i, jm]
            (isfinite(u_up) && isfinite(u_dn)) || continue
            worst = max(worst, abs(u_up - conj(u_dn)))
            ncheck += 1
        end
        @test ncheck ≥ 5
        @test worst ≤ 1e-6
    end

    @testset "FF.2.3: step-(iii) BVP fill reproduces the spine from exact anchors" begin
        # The step-(iii) primitive is a `bvp_solve` over a smooth run
        # bridged between known anchors.  Given EXACT anchors — the
        # spine values at the ends of an imaginary-axis sub-segment —
        # the fill must reproduce the spine in the interior.
        za, zb = 2.0im, 8.0im
        fill_sol = bvp_solve(bvp_f, bvp_∂f_∂u, za, zb,
                             spine(za)[1], spine(zb)[1];
                             N = 60, initial_guess = leading, maxiter = 25)
        worst = 0.0
        for y in (3.0, 4.0, 5.0, 6.0, 7.0)
            worst = max(worst, abs(fill_sol(im * y)[1] - spine(im * y)[1]))
        end
        @test worst ≤ 1e-8
    end

end # @testset FW 2011 Fig 4.1 steps (ii)+(iii)

# ----------------------------------------------------------------------
# Mutation-proof procedure (steps (ii)+(iii) testsets).
#
#   M1 — break the step-(ii)/step-(i) link.  In FF.2.1, compare
#        `egs0.u_grid` against `spine(im*ys[j] + 1.0)` (shift the spine
#        evaluation point by 1).  `worst` jumps far above 1e-6 — the
#        test genuinely pins that the two methods describe the *same*
#        solution at the *same* points.
#   M2 — break conjugacy.  In FF.2.2, compare against `egs0.u_grid[i,j]`
#        (the same cell) instead of the reflected `[i,jm]`; `worst`
#        becomes a non-trivial cell-to-cell difference and the ≤1e-6
#        bound fails on the off-axis cells.
#   M3 — wrong anchors.  In FF.2.3, pass `leading(za), leading(zb)`
#        (the crude asymptote) instead of the exact `spine` values;
#        the short segment does not smooth the BC error out and `worst`
#        rises well above 1e-8.
# Each verified RED then reverted.
# ----------------------------------------------------------------------
