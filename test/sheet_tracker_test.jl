# test/sheet_tracker_test.jl — Phase 14 / bead `padetaylor-grc`.
#
# SheetTracker: PVI ζ-plane RHS + winding-number primitives per FFW
# 2017 §2.2 (`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:144`
# the transformed PVI eq.; :163-189 circumambulation + sheet indexing).
#
# Acceptance: (a) RHS hand-pinned at a sample point with degenerate
# parameters (α=β=γ=δ=0); (b) IC compatibility with the PV transform
# (same map z↔ζ); (c) end-to-end direct-vs-transformed agreement at one
# Padé step; (d) winding primitives on closed loops + open paths;
# (e) sheet-index conversion.

using Test
using PadeTaylor

@testset "SheetTracker (Phase 14): PVI ζ-plane + winding primitives" begin

    # -----------------------------------------------------------------
    # PVI transformed RHS — hand-pinned at α=β=γ=δ=0 (degenerate)
    # -----------------------------------------------------------------
    @testset "ST.1.1: PVI ζ-plane RHS at (ζ=log 2, w=1/2, w'=1), α=β=γ=δ=0" begin
        # Hand-derivation at degenerate parameters (α=β=γ=δ=0 ⇒ third
        # term vanishes):
        #   1/w + 1/(w-1) + 1/(w-e^ζ) = 2 - 2 - 2/3 = -2/3
        #   first  = (1/2)(-2/3)(1)² = -1/3
        #   e^ζ/(e^ζ-1) + e^ζ/(w-e^ζ) = 2 - 4/3 = 2/3
        #   second = -(2/3)(1) = -2/3
        #   total  = -1/3 - 2/3 = -1
        # See module docstring for the equation; FFW2017...md:144.
        rhs = pVI_transformed_rhs(0.0, 0.0, 0.0, 0.0)
        ζ = log(2.0) + 0.0im
        w = 0.5 + 0.0im
        wp = 1.0 + 0.0im
        w_actual = rhs(ζ, w, wp)
        @test isapprox(real(w_actual), -1.0; atol = 1e-14)
        @test isapprox(imag(w_actual),  0.0; atol = 1e-14)
    end

    @testset "ST.1.2: PVI ζ-plane RHS at (ζ=log 2, w=1/2, w'=1), α=β=γ=δ=1" begin
        # Non-degenerate parameters.  Hand-derivation of the third term:
        #   w(w-1)(w-e^ζ) = (1/2)(-1/2)(-3/2) = 3/8
        #   (e^ζ-1)² = 1
        #   α + β e^ζ/w² + γ(e^ζ-1)/(w-1)² + δ e^ζ(e^ζ-1)/(w-e^ζ)²
        #     = 1 + 1·2/(1/4) + 1·1/(1/4) + 1·2·1/(9/4)
        #     = 1 + 8 + 4 + 8/9
        #     = 13 + 8/9 = 125/9
        #   third = (3/8)(125/9) = 125/24
        #   total = -1 + 125/24 = 101/24
        rhs = pVI_transformed_rhs(1.0, 1.0, 1.0, 1.0)
        w_actual = rhs(log(2.0) + 0.0im, 0.5 + 0.0im, 1.0 + 0.0im)
        @test isapprox(real(w_actual), 101.0 / 24.0; atol = 1e-14)
        @test isapprox(imag(w_actual),         0.0; atol = 1e-14)
    end

    # -----------------------------------------------------------------
    # End-to-end PVI direct-vs-transformed agreement (one Padé step)
    # -----------------------------------------------------------------
    @testset "ST.1.3: PVI end-to-end direct-vs-transformed agreement" begin
        # PVI parameters: pick a set well away from singularities for
        # the chosen step.  Stay clear of w ∈ {0, 1, e^ζ}.
        α, β, γ, δ = 0.5, -0.3, 0.2, -0.1

        # Direct PVI RHS (FFW2017...md:35):
        #   u'' = (1/2)(1/u + 1/(u-1) + 1/(u-z))(u')²
        #         - (1/z + 1/(z-1) + 1/(u-z))(u')
        #         + u(u-1)(u-z)/(z²(z-1)²) · (α + β z/u² + γ(z-1)/(u-1)² + δ z(z-1)/(u-z)²)
        f_direct = (z, u, up) -> begin
            zm1 = z - 1
            um1 = u - 1
            umz = u - z
            first  = (1/u + 1/um1 + 1/umz) * up^2 / 2
            second = -(1/z + 1/zm1 + 1/umz) * up
            param  = α + β*z/u^2 + γ*zm1/um1^2 + δ*z*zm1/umz^2
            third  = u * um1 * umz / (z * zm1)^2 * param
            return first + second + third
        end

        z₀  = 2.0 + 0.0im        # away from {0, 1}
        u₀  = 0.4 + 0.0im        # away from {0, 1, z}
        up₀ = -0.1 + 0.0im

        # Step toward z = 2.05 (small step; still away from singularities).
        z_target = 2.05 + 0.0im
        ζ_target = log(z_target)

        prob_direct = PadeTaylorProblem(f_direct, (u₀, up₀),
                                         (z₀, z_target); order = 30)
        sol_direct = path_network_solve(prob_direct,
                                         ComplexF64[z_target];
                                         h = 0.5, max_steps_per_target = 50)
        u_direct  = sol_direct.grid_u[1]
        up_direct = sol_direct.grid_up[1]

        # Transformed: PVI in ζ-plane.  Reuses PV's z↔ζ map (identical:
        # u(z) = w(ζ), z = e^ζ, so w = u and w' = z u').
        f_trans = pVI_transformed_rhs(α, β, γ, δ)
        ζ₀, w₀_, wp₀_ = pV_z_to_ζ(z₀, u₀, up₀)
        prob_trans = PadeTaylorProblem(f_trans, (w₀_, wp₀_),
                                        (ζ₀, ζ_target); order = 30)
        sol_trans = path_network_solve(prob_trans,
                                        ComplexF64[ζ_target];
                                        h = 0.5, max_steps_per_target = 50)
        w_trans  = sol_trans.grid_u[1]
        wp_trans = sol_trans.grid_up[1]

        # Invert to z-plane via PV's inverse (same map).
        _, u_rec, up_rec = pV_ζ_to_z(ζ_target, w_trans, wp_trans)

        @test isfinite(u_direct) && isfinite(u_rec)
        @test abs(u_direct  - u_rec)  ≤ 1e-10
        @test abs(up_direct - up_rec) ≤ 1e-9
    end

    # -----------------------------------------------------------------
    # Winding-number primitives
    # -----------------------------------------------------------------
    @testset "ST.1.4: winding_delta — sign + normalisation" begin
        # Step from +1 to +i around origin: Δθ = π/2 (counterclockwise).
        @test isapprox(winding_delta(1.0 + 0.0im, 0.0 + 1.0im, 0.0im),
                       π/2; atol = 1e-14)
        # Reverse: Δθ = -π/2.
        @test isapprox(winding_delta(0.0 + 1.0im, 1.0 + 0.0im, 0.0im),
                       -π/2; atol = 1e-14)
        # Branch-cut normalisation: step from just-above to just-below
        # the negative real axis.  Julia's principal `angle` is
        # discontinuous at θ = π (jumping from +π to -π); the raw
        # difference is ≈ -2π.  The normalisation to (-π, π] maps this
        # to a small SAME-SIDE wrap (≈ +0.02 here), under the implicit
        # assumption that path steps never enclose the branch point in
        # a single step.  This is the contract: callers must walk with
        # step << distance-to-branch.  Steps that ACTUALLY enclose
        # produce one full revolution of accounting error.
        z_old =  -1.0 + 0.01im
        z_new =  -1.0 - 0.01im
        Δθ = winding_delta(z_old, z_new, 0.0im)
        @test abs(Δθ) < 0.05      # normalised to a small value (~0.02)
    end

    @testset "ST.1.5: accumulate_winding — closed loops" begin
        # Counterclockwise square around origin: +2π winding.
        ccw = ComplexF64[1.0, 0.0 + 1.0im, -1.0, 0.0 - 1.0im, 1.0]
        winding_ccw = accumulate_winding(ccw, 0.0im)
        @test length(winding_ccw) == 5
        @test winding_ccw[1] == 0.0
        @test isapprox(winding_ccw[end], 2π; atol = 1e-12)

        # Clockwise: -2π.
        cw = ComplexF64[1.0, 0.0 - 1.0im, -1.0, 0.0 + 1.0im, 1.0]
        winding_cw = accumulate_winding(cw, 0.0im)
        @test isapprox(winding_cw[end], -2π; atol = 1e-12)

        # Non-enclosing path (small loop entirely in upper half-plane).
        non_enc = ComplexF64[2.0 + 1.0im, 3.0 + 1.0im,
                              3.0 + 2.0im, 2.0 + 2.0im,
                              2.0 + 1.0im]
        winding_ne = accumulate_winding(non_enc, 0.0im)
        @test abs(winding_ne[end]) ≤ 1e-12   # ≈ 0
    end

    @testset "ST.1.6: sheet_index conversion" begin
        # +2π → +1.  -2π → -1.  +4π → +2.  Sub-2π → 0 (rounds).
        @test sheet_index( 2π)         == +1
        @test sheet_index(-2π)         == -1
        @test sheet_index( 4π)         == +2
        @test sheet_index(  0.0)       ==  0
        @test sheet_index(  π / 2)     ==  0   # < π → rounds to 0
        @test sheet_index(  3π / 2)    == +1   # ≥ π → rounds up

        # Composed with accumulate_winding: ccw loop around origin
        # produces sheet +1, clockwise -1.
        ccw = ComplexF64[1.0, 0.0 + 1.0im, -1.0, 0.0 - 1.0im, 1.0]
        @test sheet_index(accumulate_winding(ccw, 0.0im)[end]) == +1
    end

    # -----------------------------------------------------------------
    # PVI ζ-plane branch-point lattice — structural check
    # -----------------------------------------------------------------
    @testset "ST.1.7: branch points at ζ = 2π·i·k blow up RHS magnitude" begin
        # FFW2017...md:141: the ζ = 2π·i·k branch lattice are fixed
        # singularities of eq. (3) where e^ζ = 1 ⇒ (e^ζ - 1) = 0.  In
        # exact arithmetic the RHS is infinite; in Float64, `exp(2π·im)`
        # is `1 + O(eps)·im` (the sin component is ≈ -2.45e-16 from
        # the rounded `2π`), so the denominator is O(eps) and the RHS
        # magnitude blows up to ~1e30 but stays technically finite.
        # The downstream Padé-Taylor stepper will see this in the
        # Taylor coefficients and throw — fail-loud per CLAUDE.md Rule 1.
        # We assert magnitude here (not non-finite) because the
        # floating-point sentinel for "branch point" is "huge", not "inf".
        rhs = pVI_transformed_rhs(1.0, 1.0, 1.0, 1.0)
        ζ_branch = 2π * im
        w = 0.5 + 0.0im
        wp = 1.0 + 0.0im
        w_at_branch = rhs(ζ_branch, w, wp)
        @test abs(w_at_branch) > 1e10    # blown up

        # Near (but not on) the branch: finite + bounded.
        ζ_near = 2π * im + 0.1 + 0.0im
        w_near = rhs(ζ_near, w, wp)
        @test isfinite(real(w_near)) && isfinite(imag(w_near))
        @test abs(w_near) < 100          # not blown up
    end

end # @testset SheetTracker

# Mutation-proof procedure (verified before commit, 2026-05-13):
#
#   Mutation O  --  in `pVI_transformed_rhs`, swap the `1/2` factor on
#     `first` (the `(dw/dζ)²` term) to `1/3`.  Verified bite: ST.1.1
#     line 35 RED (hand-pin -1 fails), ST.1.2 line 51 RED (101/24
#     fails), ST.1.3 lines 110-111 RED (end-to-end disagreement).  The
#     RHS hand-pin tests + the end-to-end transform-vs-direct check
#     each independently catch the algebraic error.
#
#   Mutation P  --  in `winding_delta`, drop the wrap-to-(-π, π]
#     normalisation entirely (return raw `Δθ`).  Verified bite: ST.1.4
#     line 136 (near-branch-cut hop's normalised value blows up to
#     ≈-2π absolute); ST.1.5 lines 145 + 150 (closed-loop CCW + CW
#     accumulations both wrong by ±2π due to the unnormalized
#     intermediate `Δθ`s near the cut); ST.1.6 line 172 (downstream
#     `sheet_index` derives from accumulate_winding).
#
#   Mutation Q  --  in `sheet_index`, replace `round` with `floor`.
#     Verified bite: ST.1.6 line 167 RED — sheet_index(3π/2) reads 0
#     instead of +1 (floor(0.75) = 0).  Other ST.1.6 assertions pass
#     because their inputs are at integer multiples of 2π where round
#     and floor agree.  One assertion catching the mutation is enough
#     to confirm load-bearing; the test could be strengthened with
#     additional sub-2π inputs but the principle is established.
#
# Restoration: all three mutations restored before commit; GREEN at 1311.
