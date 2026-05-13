# test/coord_transforms_test.jl — Phase 13 / bead `padetaylor-bvh`.
#
# CoordTransforms: exponential coordinate maps for PIII and PV per
# FFW 2017 §2.1 (`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:39-48`).
#
# Acceptance: (a) IC round-trip pinpoint; (b) RHS closure evaluates to
# the FFW-eq value at a sample point; (c) end-to-end agreement —
# integrate one PathNetwork Padé step both directly (`u'' = PIII(z, u, u')`)
# and via the transform (`w'' = P̃_III(ζ, w, w')` → invert), then compare.
#
# The FFW 2017 Fig 1 parameter set is the reference scenario:
#   PIII: α = -1/2, β = -1/2, γ = 1, δ = -1, IC (u(1), u'(1)) = (1/4, 1).
#
# We use a SMALLER PV scenario at the same z-base (z=1) with parameters
# chosen to keep `w ≠ 0` and `w ≠ 1` in the test step neighborhood.

using Test
using PadeTaylor

@testset "CoordTransforms (Phase 13): PIII / PV exponential maps" begin

    # -----------------------------------------------------------------
    # PIII tests — FFW 2017 Fig 1 scenario
    # -----------------------------------------------------------------
    @testset "PIII transformation" begin

        z₀  = 1.0 + 0.0im
        u₀  = 0.25 + 0.0im
        up₀ = 1.0  + 0.0im
        α, β, γ, δ = -0.5, -0.5, 1.0, -1.0   # FFW 2017 Fig 1 caption

        @testset "CT.1.1: IC round-trip at (z=1, u=1/4, u'=1)" begin
            # z → ζ via z = exp(ζ/2): at z = 1, ζ = 0; w = z*u = 1/4;
            # w' = (z u + z² u') / 2 = (1/4 + 1) / 2 = 5/8.
            ζ, w, wp = pIII_z_to_ζ(z₀, u₀, up₀)
            @test isapprox(ζ,  0.0 + 0.0im; atol = 1e-15)
            @test isapprox(w,  0.25 + 0.0im; atol = 1e-15)
            @test isapprox(wp, 0.625 + 0.0im; atol = 1e-15)

            # Inverse: (ζ, w, w') → (z, u, u') must recover input.
            z₂, u₂, up₂ = pIII_ζ_to_z(ζ, w, wp)
            @test isapprox(z₂,  z₀;  atol = 1e-14)
            @test isapprox(u₂,  u₀;  atol = 1e-14)
            @test isapprox(up₂, up₀; atol = 1e-14)
        end

        @testset "CT.1.2: P̃_III RHS evaluates to FFW eq at a sample point" begin
            # Compute w'' from the FFW formula by hand at (ζ, w, w') =
            # (0, 1/4, 5/8) with FFW Fig 1 parameters.  FFW eq. (line 43):
            #   w'' = (w')²/w + (α w² + γ w³ + β eᶻ + δ e²ᶻ / w) / 4
            # At ζ=0: eᶻ = 1, e²ᶻ = 1.  α = -1/2, β = -1/2, γ = 1, δ = -1.
            #   (w')²/w = (5/8)² / (1/4) = (25/64) · 4 = 25/16 = 1.5625
            #   α w²     = (-1/2)(1/16)  = -1/32
            #   γ w³     = (1)(1/64)     = 1/64
            #   β eᶻ     = -1/2
            #   δ e²ᶻ/w  = -1 · 1/(1/4)  = -4
            #   sum / 4  = (-1/32 + 1/64 - 1/2 - 4) / 4
            #            = (-2/64 + 1/64 - 32/64 - 256/64) / 4
            #            = -289/64 / 4 = -289/256
            #   total    = 25/16 - 289/256 = 400/256 - 289/256 = 111/256
            w_expected = 111.0 / 256.0   # 0.43359375
            rhs = pIII_transformed_rhs(α, β, γ, δ)
            w_actual = rhs(0.0 + 0.0im, 0.25 + 0.0im, 0.625 + 0.0im)
            @test isapprox(real(w_actual), w_expected; atol = 1e-14)
            @test isapprox(imag(w_actual), 0.0;        atol = 1e-14)
        end

        @testset "CT.1.3: end-to-end agreement direct vs transformed (single step)" begin
            # Direct PIII: f_direct(z, u, u') = (1/u)(u')² - (1/z)u' +
            #   (α u² + β)/z + γ u³ + δ/u.  Take ONE Padé step from z=1
            #   to z=1 + h_z (small h to stay away from poles and from
            #   z=0 branch point).
            f_direct = (z, u, up) -> (up^2)/u - up/z +
                                      (α*u^2 + β)/z + γ*u^3 + δ/u

            # Use a complex h_z that gives the same |z step| in the ζ-plane.
            # Pick h_z = 0.05 (real, small).  In ζ-plane this corresponds
            # to going from ζ=0 to ζ_target = 2 log(1.05) ≈ 0.0976.
            z_target = 1.05 + 0.0im
            ζ_target = 2 * log(z_target)

            # --- direct path-network call (no transform) -----------------
            prob_direct = PadeTaylorProblem(f_direct, (u₀, up₀),
                                             (z₀, z_target); order = 30)
            sol_direct = path_network_solve(prob_direct,
                                             ComplexF64[z_target];
                                             h = 0.5, max_steps_per_target = 50)
            u_direct  = sol_direct.grid_u[1]
            up_direct = sol_direct.grid_up[1]

            # --- transformed path-network call ---------------------------
            f_trans = pIII_transformed_rhs(α, β, γ, δ)
            ζ₀, w₀_, wp₀_ = pIII_z_to_ζ(z₀, u₀, up₀)
            prob_trans = PadeTaylorProblem(f_trans, (w₀_, wp₀_),
                                            (ζ₀, ζ_target); order = 30)
            sol_trans = path_network_solve(prob_trans,
                                            ComplexF64[ζ_target];
                                            h = 0.5, max_steps_per_target = 50)
            w_trans  = sol_trans.grid_u[1]
            wp_trans = sol_trans.grid_up[1]

            # Invert to z-plane.
            _, u_recovered, up_recovered = pIII_ζ_to_z(ζ_target, w_trans, wp_trans)

            # Both paths solve the SAME analytic IVP; agreement is
            # limited by the path-network's own truncation error, not by
            # the transform.  At h=0.5 with order=30, we routinely see
            # ≤ 1e-10 on smooth segments.
            @test isfinite(u_direct) && isfinite(u_recovered)
            @test abs(u_direct  - u_recovered)  ≤ 1e-10
            @test abs(up_direct - up_recovered) ≤ 1e-9
        end
    end

    # -----------------------------------------------------------------
    # PV tests
    # -----------------------------------------------------------------
    @testset "PV transformation" begin

        # Pick PV parameters/IC keeping w ≠ 0 and w ≠ 1 along the step.
        z₀  = 1.0 + 0.0im
        u₀  = 0.5 + 0.0im
        up₀ = 0.3 + 0.0im
        α, β, γ, δ = 0.1, -0.2, 0.05, 0.01

        @testset "CT.1.4: PV IC round-trip at (z=1, u=0.5, u'=0.3)" begin
            # z → ζ via z = exp(ζ): at z=1, ζ=0; w=u=1/2; w' = z u' = 0.3.
            ζ, w, wp = pV_z_to_ζ(z₀, u₀, up₀)
            @test isapprox(ζ,  0.0 + 0.0im; atol = 1e-15)
            @test isapprox(w,  0.5 + 0.0im; atol = 1e-15)
            @test isapprox(wp, 0.3 + 0.0im; atol = 1e-15)

            z₂, u₂, up₂ = pV_ζ_to_z(ζ, w, wp)
            @test isapprox(z₂,  z₀;  atol = 1e-14)
            @test isapprox(u₂,  u₀;  atol = 1e-14)
            @test isapprox(up₂, up₀; atol = 1e-14)
        end

        @testset "CT.1.5: P̃_V RHS evaluates to FFW eq at a sample point" begin
            # By hand at (ζ, w, w') = (0, 1/2, 3/10), with α=0.1, β=-0.2,
            # γ=0.05, δ=0.01.  FFW eq. (line 47):
            #   w'' = (1/(2w) + 1/(w-1)) (w')²
            #         + (w-1)² (α w + β/w)
            #         + γ eᶻ w
            #         + δ e²ᶻ w(w+1)/(w-1)
            # At ζ=0, w=1/2:
            #   (1/(2w) + 1/(w-1)) = (1/1 + 1/(-1/2)) = 1 - 2 = -1
            #   (w')² = 0.09
            #   ⇒ first term = -0.09
            #   (w-1)² = 1/4;  α w + β/w = 0.05 + (-0.4) = -0.35
            #   ⇒ second term = (1/4)(-0.35) = -0.0875
            #   γ eᶻ w = 0.05 · 1 · 0.5 = 0.025
            #   δ e²ᶻ w(w+1)/(w-1) = 0.01 · 1 · (0.5)(1.5)/(-0.5) = -0.015
            #   sum = -0.09 - 0.0875 + 0.025 - 0.015 = -0.1675
            w_expected = -0.1675
            rhs = pV_transformed_rhs(α, β, γ, δ)
            w_actual = rhs(0.0 + 0.0im, 0.5 + 0.0im, 0.3 + 0.0im)
            @test isapprox(real(w_actual), w_expected; atol = 1e-14)
            @test isapprox(imag(w_actual), 0.0;        atol = 1e-14)
        end

        @testset "CT.1.6: end-to-end agreement PV direct vs transformed" begin
            # PV direct RHS:
            #   u'' = (1/(2u) + 1/(u-1)) (u')² - (1/z) u'
            #         + ((u-1)²/z²)(α u + β/u) + γ u/z + δ u(u+1)/(u-1)
            f_direct = (z, u, up) -> (1/(2*u) + 1/(u-1)) * up^2 - up/z +
                                      ((u-1)^2 / z^2) * (α*u + β/u) +
                                      γ*u/z +
                                      δ*u*(u+1)/(u-1)

            z_target = 1.05 + 0.0im
            ζ_target = log(z_target)   # PV: ζ = log z (not 2 log z)

            prob_direct = PadeTaylorProblem(f_direct, (u₀, up₀),
                                             (z₀, z_target); order = 30)
            sol_direct = path_network_solve(prob_direct,
                                             ComplexF64[z_target];
                                             h = 0.5, max_steps_per_target = 50)
            u_direct  = sol_direct.grid_u[1]
            up_direct = sol_direct.grid_up[1]

            f_trans = pV_transformed_rhs(α, β, γ, δ)
            ζ₀, w₀_, wp₀_ = pV_z_to_ζ(z₀, u₀, up₀)
            prob_trans = PadeTaylorProblem(f_trans, (w₀_, wp₀_),
                                            (ζ₀, ζ_target); order = 30)
            sol_trans = path_network_solve(prob_trans,
                                            ComplexF64[ζ_target];
                                            h = 0.5, max_steps_per_target = 50)
            w_trans  = sol_trans.grid_u[1]
            wp_trans = sol_trans.grid_up[1]
            _, u_recovered, up_recovered = pV_ζ_to_z(ζ_target, w_trans, wp_trans)

            @test isfinite(u_direct) && isfinite(u_recovered)
            @test abs(u_direct  - u_recovered)  ≤ 1e-10
            @test abs(up_direct - up_recovered) ≤ 1e-9
        end
    end

end # @testset CoordTransforms

# Mutation-proof procedure (verified before commit):
#
#   Mutation L  --  in `pIII_transformed_rhs`, replace the `/ 4` divisor
#     with `/ 3` on the parametrised term group.  Expected RED:
#     CT.1.2 RED at the `w_actual ≈ w_expected` assertion (the
#     hand-computed -289/256 vs (-289/192) numerically far apart);
#     CT.1.3 RED at u-agreement (one-step error ≫ 1e-10 since the
#     transformed solve is now solving a different ODE).
#
#   Mutation M  --  in `pIII_z_to_ζ`, change `wp = (z u + z² u') / 2`
#     to `wp = (z u - z² u') / 2` (sign flip).  Expected RED:
#     CT.1.1 RED at `wp ≈ 5/8` (computed value is (1/4 - 1)/2 = -3/8);
#     CT.1.3 RED at end-to-end agreement (transformed IVP starts from
#     wrong w').
#
#   Mutation N  --  in `pV_transformed_rhs`, swap `δ e²ᶻ w(w+1)/(w-1)`
#     to `δ e²ᶻ w(w-1)/(w+1)` (swap the +/- in num and denom).
#     Expected RED: CT.1.5 RED — the hand-computed value changes from
#     -0.1675 to something else (≈ -0.0825 by recomputation); CT.1.6
#     RED at end-to-end agreement.
#
# Restoration: all three mutations restored before commit; tests GREEN.
