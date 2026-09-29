# test/ivp_bvp_hybrid_branch_test.jl
#
# Bead padetaylor-w80i — the branch of z^{1/3} on the upper part of the
# FFW 2017 Fig 5 pole-free sector.
#
# FFW md:222 (references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
# FFW2017_painleve_riemann_surfaces_preprint.md:222) puts the tronquée P_III
# sector at -3π/4 < arg z < 9π/4, i.e. -3π/2 < Im ζ < 9π/2 under z = exp(ζ/2).
# The sector is 3π wide, so for Im ζ > 2π the z-frame point wraps onto the
# principal argument range.  Before the fix, `solve_pole_free_hybrid` formed
# z = exp(ζ/2) and called the IC series on the principal root: at FFW's own
# upper IC point z₁ = 30e^{13πi/6} (md:243) that gives |u − u_FFW| = 5.39
# (worklog 087), a different sheet.  The fix threads the continuous
# argument Im ζ / 2 through as `argz`.
#
# Why the ODE residual does NOT discriminate the branches (checked in BR.3):
# the ansatz u = s + a_1 s^{-1} + … is a Laurent series in s = z^{1/3}
# itself, and z = s³ on every branch, so a formal solution in s is a formal
# solution for every choice of cube root.  The three branches are three
# distinct tronquée solutions (FFW md:222: "any of the branches of z^{1/3}").
# Only FFW's published value at z₁ (md:243) — i.e. analytic continuation —
# selects the right one.  So the load-bearing oracle here is md:243.

using Test
using PadeTaylor
using PadeTaylor.IVPBVPHybrid: _call_asymptotic_on_branch, _branch_arg

@testset "IVPBVPHybrid — z^{1/3} branch on the upper sector (bead w80i)" begin
    β, δ = -1/20, -1.0
    # FFW md:243 verbatim.
    u_z1_ffw  = -2.000735432319 + 2.376177147900im
    up_z1_ffw = -5.939523100e-3 + 3.402038641e-2im
    u_z2_ffw  =  2.384379236170 - 1.993845650158im
    ζ1 = complex(2log(30.0), 13π/3)     # z₁ = 30e^{13πi/6}; Im ζ > 2π
    ζ2 = complex(2log(30.0), -4π/3)     # z₂ = 30e^{-2πi/3}; principal
    z1 = exp(ζ1 / 2);  z2 = exp(ζ2 / 2)

    @testset "BR.1: argz selects the continued branch; matches FFW u(z₁)" begin
        u, up = pIII_asymptotic_ic(z1; argz = 13π/6, n_terms = 2, β = β, δ = δ)
        # 8.4e-6 measured: the dropped a_3 z^{-7/3} term (docstring).
        @test abs(u  - u_z1_ffw)  < 5e-5
        @test abs(up - up_z1_ffw) < 5e-6
        # The principal root at the same z is a different sheet (5.39 away).
        u_pr, _ = pIII_asymptotic_ic(z1; n_terms = 2, β = β, δ = δ)
        @test abs(u_pr - u_z1_ffw) > 5.0
    end

    @testset "BR.2: argz validation (fail loud)" begin
        # argz must agree with z mod 2π.
        @test_throws ArgumentError pIII_asymptotic_ic(z1; argz = 13π/6 - 2π/3)
        # argz must lie in the FFW sector (-3π/4, 9π/4).
        @test_throws ArgumentError pIII_asymptotic_ic(30.0 * cis(9π/4 + 0.1);
                                                      argz = 9π/4 + 0.1)
        # On the principal slice argz = angle(z) reproduces the default.
        u_d, up_d = pIII_asymptotic_ic(z2; β = β, δ = δ)
        u_a, up_a = pIII_asymptotic_ic(z2; argz = angle(z2), β = β, δ = δ)
        @test isapprox(u_a,  u_d;  rtol = 1e-14)
        @test isapprox(up_a, up_d; rtol = 1e-14)
        @test abs(u_d - u_z2_ffw) < 5e-5
    end

    @testset "BR.3: ODE residual is small on BOTH branches (not an oracle)" begin
        # u = s + a1/s (a_2 = 0 at δ = −1), s = z^{1/3}; closed-form
        # derivatives in s via d/dz = (1/(3s²)) d/ds.
        a1 = -β / 3
        function pIII_residual(s)
            z   = s^3
            u   = s + a1 / s
            up  = (s^-2 - a1 * s^-4) / 3
            upp = (-2 * s^-3 + 4 * a1 * s^-5) / (9 * s^2)
            return upp - (up^2 / u - up / z + (u^2 + β) / z + δ / u)
        end
        s_cont = cbrt(30.0) * cis(13π / 18)   # continued branch at z₁
        s_pr   = z1^(1/3)                      # principal branch at z₁
        @test abs(pIII_residual(s_cont)) < 1e-5
        @test abs(pIII_residual(s_pr))   < 1e-5
        @test abs((s_pr + a1 / s_pr) - (s_cont + a1 / s_cont)) > 5.0
    end

    @testset "BR.4: driver helper passes Im ζ/2 as argz" begin
        pp = PainleveProblem(:III; α = 1.0, β = β, γ = 0.0, δ = δ,
                             u0 = u_z2_ffw, up0 = 0.0im,
                             zspan = (z2, z2 + 1.0), order = 30)
        @test _branch_arg(pp, ζ1) == 13π/6
        fn_kw = (z; argz = nothing) ->
            pIII_asymptotic_ic(z; argz, n_terms = 2, β = β, δ = δ)
        u, _ = _call_asymptotic_on_branch(fn_kw, z1, _branch_arg(pp, ζ1), "top")
        @test abs(u - u_z1_ffw) < 5e-5
        # A callable that cannot take `argz` must not silently get the
        # principal branch off the principal slice.
        fn_plain = z -> pIII_asymptotic_ic(z; n_terms = 2, β = β, δ = δ)
        @test_throws ArgumentError _call_asymptotic_on_branch(
            fn_plain, z1, _branch_arg(pp, ζ1), "top")
        # On the principal slice the plain callable is still accepted.
        u2, _ = _call_asymptotic_on_branch(fn_plain, z2, _branch_arg(pp, ζ2), "bot")
        @test abs(u2 - u_z2_ffw) < 5e-5
    end

    @testset "BR.5: full driver on FFW's upper sector Im ζ up to 13π/3" begin
        # FFW md:222 sector, top edge placed so ζ_top sits at FFW's z₁.
        up_z2_ffw = 6.050817704e-3 + 3.398020750e-2im
        pp = PainleveProblem(:III; α = 1.0, β = β, γ = 0.0, δ = δ,
                             u0 = u_z2_ffw, up0 = up_z2_ffw,
                             zspan = (z2, z2 + 1.0), order = 30)
        sector = (im_lo = -3π/2 + 0.05, im_hi = 13π/3 + 0.01,
                  re_anchor = 2log(30.0), re_extent = 1.0)
        kw = (pfs_kwargs = (; h = 0.4, step_size_policy = :adaptive_ffw,
                            adaptive_tol = 1e-10, k_conservative = 1e-3,
                            max_rescales = 50, max_steps_per_target = 500),
              bvp_kwargs = (; N = 30, tol = 1e-10, max_iter = 30),
              n_slices = 5, glue_tol = 1e-8)
        fn_kw = (z; argz = nothing) ->
            pIII_asymptotic_ic(z; argz, n_terms = 2, β = β, δ = δ)
        fn_plain = z -> pIII_asymptotic_ic(z; n_terms = 2, β = β, δ = δ)
        @test_throws ArgumentError solve_pole_free_hybrid(pp, sector, fn_plain; kw...)
        sol = solve_pole_free_hybrid(pp, sector, fn_kw; kw...)
        # Reference: the continued-branch series in ζ, u ~ e^{ζ/6}(1 + a_1
        # e^{-ζ/3}) (FFW md:222), w = e^{ζ/2} u.  Measured relative gap
        # 2.2e-3 at every height (v1 hybrid accuracy, not the branch);
        # a principal-branch IC gives 1.74 ≈ |1 − e^{2πi/3}| = √3 near
        # the top edge (worklog 087).
        a1 = -β / 3
        wref(ζ) = exp(ζ / 2) * (exp(ζ / 6) + a1 * exp(-ζ / 6))
        for y in (0.0, 8.0, 13.0)
            ζ = complex(2log(30.0) - 0.3, y)
            @test abs(sol(ζ)[1] - wref(ζ)) / abs(wref(ζ)) < 1e-2
        end
    end
end
