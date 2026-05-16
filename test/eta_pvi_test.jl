# test/eta_pvi_test.jl — Phase 14 follow-on / bead `padetaylor-riu`.
#
# Step A3 of the 11-step FFW 2017 reproduction plan: η-plane PVI RHS
# (FFW eq. 5) and the `z ↔ η` IC helpers, plus a `:transformed_eta`
# frame for `PainleveProblem(:VI; ...)`.
#
# Ground truth: `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
# FFW2017_painleve_riemann_surfaces_preprint.md:146-162`.  Specifically
# md:154 is FFW eq. (5), the η-plane PVI equation; md:157 records the
# branch-point-free region `Re η < log(2π)`; md:195 gives FFW Fig 2's
# IC point: `(α,β,γ,δ) = (4,-4,8,-8)`, `u(10) = 0.429534600325223`,
# `u'(10) = -1.61713114374804e-3`.
#
# Acceptance (mirrors `sheet_tracker_test.jl`'s pattern):
#
#   ET.1.1 — RHS hand-pin at α=β=γ=δ=0 (third term vanishes, leaving
#            only first + second).  Closed-form expected value
#            derived from the algebra inline.
#   ET.1.2 — IC round-trip `pVI_z_to_η ∘ pVI_η_to_z = id` (and
#            symmetrically) to machine precision at the FFW Fig 2 IC
#            point.
#   ET.1.3 — End-to-end direct (z-plane) vs η-transformed agreement on
#            one Padé step.  Same shape as ST.1.3 and CT.1.3.
#   ET.1.4 — Branch-point-free region check (md:157): RHS finite at
#            sample points with `Re η < log(2π)`; transitions to
#            huge-but-finite at sample points outside.

using Test
using PadeTaylor

@testset "η-plane PVI (Phase 14 follow-on): RHS + IC helpers + frame wiring" begin

    # -----------------------------------------------------------------
    # ET.1.1 — RHS hand-pin at α=β=γ=δ=0
    # -----------------------------------------------------------------
    #
    # Pick `η = log(log 2)` so that the nested exponential collapses to
    # a clean integer: `E := exp(exp(η)) = exp(log 2) = 2`, with
    # `e^η = log 2`.  Pick `v = 1/2`, `v' = 1`.  With all four
    # parameters zero, the third term of FFW eq. (5) vanishes
    # identically; the equation reduces to first + second:
    #
    #   First term:
    #     (1/v + 1/(v-1) + 1/(v-E)) (v')² / 2
    #       = (1/(1/2) + 1/(1/2 - 1) + 1/(1/2 - 2)) · 1 / 2
    #       = (2 + (-2) + (-2/3)) · 1/2
    #       = (-2/3) / 2
    #       = -1/3
    #
    #   Second term:
    #     -(e^η·E/(E-1) + e^η·E/(v-E) - 1) · v'
    #       = -((log 2)·2/1 + (log 2)·2/(-3/2) - 1) · 1
    #       = -(2 log 2 - (4/3) log 2 - 1)
    #       = -((2/3) log 2 - 1)
    #       = 1 - (2/3) log 2
    #
    #   Third term: 0.
    #
    #   Total: -1/3 + 1 - (2/3) log 2 = 2/3 - (2/3) log 2
    #        = (2/3)(1 - log 2) ≈ 0.20456854629336977
    #
    # The expected value is computed inline from the *closed-form
    # algebra*, not from a recompute of the impl formula — that's what
    # makes this a real hand-pin (Mutation M1 below confirms the
    # impl's bracketed-(-1) term is load-bearing on this assertion).
    @testset "ET.1.1: η-plane RHS at (η=log log 2, v=1/2, v'=1), α=β=γ=δ=0" begin
        rhs      = pVI_eta_transformed_rhs(0.0, 0.0, 0.0, 0.0)
        η        = log(log(2.0)) + 0.0im
        v        = 0.5 + 0.0im
        vp       = 1.0 + 0.0im
        v_actual = rhs(η, v, vp)
        expected = (2/3) * (1 - log(2))           # ≈ 0.20456854629336977
        @test isapprox(real(v_actual), expected; atol = 1e-14)
        @test isapprox(imag(v_actual), 0.0;      atol = 1e-14)
    end

    # ET.1.1bis — same point with non-trivial parameters α=β=γ=δ=1,
    # extending the hand-pin to cover the third (parameter) term.
    #
    #   Third term coefficient:
    #     v(v-1)(v-E) · e^(2η) / (E-1)²
    #       = (1/2)(-1/2)(-3/2) · (log 2)² / 1
    #       = (3/8)(log 2)²
    #
    #   Parameter bracket with α=β=γ=δ=1:
    #     1 + β·E/v² + γ(E-1)/(v-1)² + δ E(E-1)/(v-E)²
    #       = 1 + 1·2/(1/4) + 1·1/(1/4) + 1·2·1/(9/4)
    #       = 1 + 8 + 4 + 8/9
    #       = 13 + 8/9 = 125/9
    #
    #   Third term: (3/8)(log 2)² · (125/9) = (125/24)(log 2)²
    #
    #   Grand total: 2/3 - (2/3) log 2 + (125/24)(log 2)²
    #              ≈ 2.706927993784002
    @testset "ET.1.1bis: η-plane RHS at same point, α=β=γ=δ=1" begin
        rhs      = pVI_eta_transformed_rhs(1.0, 1.0, 1.0, 1.0)
        η        = log(log(2.0)) + 0.0im
        v        = 0.5 + 0.0im
        vp       = 1.0 + 0.0im
        v_actual = rhs(η, v, vp)
        expected = (2/3) * (1 - log(2)) + (125/24) * log(2)^2
        @test isapprox(real(v_actual), expected; atol = 1e-14)
        @test isapprox(imag(v_actual), 0.0;      atol = 1e-14)
    end

    # -----------------------------------------------------------------
    # ET.1.2 — IC round-trip at FFW Fig 2's IC point (md:195)
    # -----------------------------------------------------------------
    @testset "ET.1.2: pVI_z_to_η ↔ pVI_η_to_z round-trip at FFW Fig 2 IC" begin
        # FFW2017...md:195: Fig 2 caption.
        z   = 10.0 + 0.0im
        u   = 0.429534600325223 + 0.0im
        up  = -1.61713114374804e-3 + 0.0im

        η, v, vp = pVI_z_to_η(z, u, up)

        # Sanity: η should be real (z is real positive), and equal
        # log(log z) = log(log 10) ≈ 0.834.  v should be u (the
        # dependent variable is unchanged under z → η).
        @test isapprox(real(η), log(log(10.0)); atol = 1e-14)
        @test isapprox(imag(η), 0.0;            atol = 1e-14)
        @test v == u

        # Chain rule: vp = z·ζ·up = z·log(z)·up.
        @test isapprox(vp, 10.0 * log(10.0) * up; atol = 1e-14)

        # Round-trip via inverse.
        z_back, u_back, up_back = pVI_η_to_z(η, v, vp)
        @test abs(z_back  - z)  ≤ 1e-13
        @test abs(u_back  - u)  ≤ 1e-15
        @test abs(up_back - up) ≤ 1e-15

        # Symmetric round-trip (η → z → η).
        z2, u2, up2 = pVI_η_to_z(η, v, vp)
        η2, v2, vp2 = pVI_z_to_η(z2, u2, up2)
        @test abs(η2  - η)  ≤ 1e-13
        @test abs(v2  - v)  ≤ 1e-15
        @test abs(vp2 - vp) ≤ 1e-13
    end

    # -----------------------------------------------------------------
    # ET.1.3 — End-to-end direct (z-plane) vs η-transformed agreement.
    # Same Padé step from `z₀ = 5` to `z = 5.05` solved two ways; the
    # η-plane recovered `u` must match the direct `u` at the target.
    # Proves the transform is faithful: same analytic solution.
    # -----------------------------------------------------------------
    @testset "ET.1.3: PVI direct-vs-η-transformed agreement at one step" begin
        α, β, γ, δ = 0.5, -0.3, 0.2, -0.1

        # Direct PVI RHS in z-plane (FFW2017...md:35).
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

        # IC away from {0, 1}.  Pick z₀ = 5 (well inside the
        # branch-point-free η region: Re η = log log 5 ≈ 0.48 < log(2π)
        # ≈ 1.84).
        z₀  = 5.0 + 0.0im
        u₀  = 0.4 + 0.0im
        up₀ = -0.1 + 0.0im

        z_target = 5.05 + 0.0im
        η_target = log(log(z_target))      # principal branches

        prob_direct = PadeTaylorProblem(f_direct, (u₀, up₀),
                                         (z₀, z_target); order = 30)
        sol_direct  = path_network_solve(prob_direct,
                                          ComplexF64[z_target];
                                          h = 0.5,
                                          max_steps_per_target = 50)
        u_direct  = sol_direct.grid_u[1]
        up_direct = sol_direct.grid_up[1]

        # η-plane: FFW eq. (5).
        f_eta             = pVI_eta_transformed_rhs(α, β, γ, δ)
        η₀, v₀, vp₀       = pVI_z_to_η(z₀, u₀, up₀)
        prob_eta          = PadeTaylorProblem(f_eta, (v₀, vp₀),
                                              (η₀, η_target); order = 30)
        sol_eta           = path_network_solve(prob_eta,
                                                ComplexF64[η_target];
                                                h = 0.5,
                                                max_steps_per_target = 50)
        v_eta  = sol_eta.grid_u[1]
        vp_eta = sol_eta.grid_up[1]

        # Invert to z-plane.
        _, u_rec, up_rec = pVI_η_to_z(η_target, v_eta, vp_eta)

        @test isfinite(u_direct) && isfinite(u_rec)
        @test abs(u_direct  - u_rec)  ≤ 1e-10
        @test abs(up_direct - up_rec) ≤ 1e-9
    end

    # -----------------------------------------------------------------
    # ET.1.4 — Branch-point-free region check (md:157).
    # `Re η < log(2π) ≈ 1.83788` is meant to be free of η-plane branch
    # points.  Probe RHS finiteness inside, just outside, and on a
    # branch point.
    # -----------------------------------------------------------------
    @testset "ET.1.4: branch-point-free region (FFW2017...md:157)" begin
        rhs = pVI_eta_transformed_rhs(1.0, 1.0, 1.0, 1.0)
        v   = 0.5 + 0.0im
        vp  = 1.0 + 0.0im

        # Inside Re η < log(2π): finite everywhere.
        for η in (0.0 + 0.0im, 0.5 + 1.0im, 1.0 + 0.0im, 1.5 - 2.0im)
            r = rhs(η, v, vp)
            @test isfinite(real(r))
            @test isfinite(imag(r))
            @test abs(r) < 1e5   # bounded; not blown up
        end

        # Just outside Re η = log(2π) on the real axis: the k = ±1
        # branch points are at η = log(2π) + i·(±π/2), so a point at
        # `Re η = log(2π) + 0.1`, `Im η = 0` is away from them on the
        # real axis but past the boundary of the pole-free region.
        # The RHS is still finite (no branch point on the real axis
        # between |k| ≥ 1 lattice points), just larger.
        η_outside = log(2π) + 0.1 + 0.0im
        r_outside = rhs(η_outside, v, vp)
        @test isfinite(real(r_outside))
        @test isfinite(imag(r_outside))

        # On a branch point: `η = log(2π) + i·(π/2)` is the k = +1
        # lattice point (md:148, eq. 4: `η = log|2πk| + i·arg(2π·i·k)`
        # with k = +1 gives `log(2π) + i·arg(2π·i) = log(2π) + i·π/2`).
        # In exact arithmetic the RHS is ∞; in Float64 the denominator
        # `E - 1 = exp(2π·i) - 1` is `O(eps)`-small but non-zero, so
        # the RHS blows up to ~1e31 — huge but technically finite.
        # Same fail-loud-at-stepper semantics as ST.1.7.
        η_branch  = log(2π) + (π/2) * im
        r_branch  = rhs(η_branch, v, vp)
        @test abs(r_branch) > 1e20    # blown up
    end

    # -----------------------------------------------------------------
    # ET.1.5 — `:transformed_eta` frame wiring through PainleveProblem
    # -----------------------------------------------------------------
    @testset "ET.1.5: PainleveProblem(:VI; frame=:transformed_eta)" begin
        α, β, γ, δ = 4.0, -4.0, 8.0, -8.0
        z          = 10.0 + 0.0im
        u          = 0.429534600325223 + 0.0im
        up         = -1.61713114374804e-3 + 0.0im

        # Default frame: :transformed (ζ-plane).
        pp_ζ = PainleveProblem(:VI; α = α, β = β, γ = γ, δ = δ,
                               u0 = u, up0 = up, zspan = (z, 5.0 + 0.0im))
        @test pp_ζ.frame == :transformed

        # Explicit η-plane frame.
        pp_η = PainleveProblem(:VI; α = α, β = β, γ = γ, δ = δ,
                               u0 = u, up0 = up, zspan = (z, 5.0 + 0.0im),
                               frame = :transformed_eta)
        @test pp_η.frame == :transformed_eta
        @test pp_η.equation == :VI

        # The η-frame problem's IC point in solve frame must equal
        # `pVI_z_to_η(z, u, up)` (the (η, v, vp) tuple); the
        # `PadeTaylorProblem`'s IC stores `(w₀, w'₀)` = `(v, vp)`.
        η_expected, v_expected, vp_expected = pVI_z_to_η(z, u, up)
        @test pp_η.problem.zspan[1] ≈ η_expected
        @test pp_η.problem.y0[1] ≈ v_expected
        @test pp_η.problem.y0[2] ≈ vp_expected

        # The to_frame / from_frame round-trip is the identity (the
        # closures wire `pVI_z_to_η` and `pVI_η_to_z`).
        η_round, v_round, vp_round =
            pp_η.to_frame(z, u, up)
        z_back, u_back, up_back =
            pp_η.from_frame(η_round, v_round, vp_round)
        @test abs(z_back  - z)  ≤ 1e-13
        @test abs(u_back  - u)  ≤ 1e-15
        @test abs(up_back - up) ≤ 1e-15

        # Unknown frame fails loud with a helpful message.
        @test_throws ArgumentError PainleveProblem(:VI;
            α = α, β = β, γ = γ, δ = δ, u0 = u, up0 = up,
            zspan = (z, 5.0 + 0.0im), frame = :unknown_frame)
    end

end # @testset η-plane PVI

# Mutation-proof procedure (re-verified before commit, 2026-05-16,
# at worklog 041; 42 GREEN under the un-mutated impl):
#
#   Mutation M1  --  in `pVI_eta_transformed_rhs`, drop the `- 1` term
#     inside the second-term bracket (the chain-rule artefact from
#     stacking ζ = e^η on top of z = e^ζ; cf. ζ-plane eq. 3 which has
#     no such term).  Concretely: change
#       second = -(eη * E / E_m1 + eη * E / v_mE - 1) * vp
#     to
#       second = -(eη * E / E_m1 + eη * E / v_mE) * vp
#     Verified bite (4 RED of 42): ET.1.1 (1) + ET.1.1bis (1) RED on
#     the hand-pinned closed-form RHS values (the `- 1` is load-
#     bearing on the constant offset); ET.1.3 (2) RED on the end-to-
#     end direct-vs-η-transformed agreement (the η-frame RHS is no
#     longer FFW eq. 5).
#
#   Mutation M2  --  in `pVI_z_to_η`, flip the sign on the chain-rule
#     factor.  Concretely: change `vp = z * ζ * up` to `vp = -z * ζ
#     * up`.  Verified bite (6 RED of 42): ET.1.2 (3) RED on the
#     direct chain-rule check + both round-trip assertions (forward
#     mutates but inverse doesn't, so up_back = -up); ET.1.3 (2) RED
#     on the end-to-end direct-vs-η; ET.1.5 (1) RED on the
#     PainleveProblem round-trip through `to_frame`/`from_frame`.
#
#   Mutation M3  --  in `_build_VI`, dispatch `:transformed_eta` to
#     the ζ-plane RHS by mistake.  Concretely: in the
#     `:transformed_eta` branch replace `pVI_eta_transformed_rhs`
#     with `pVI_transformed_rhs` and `pVI_z_to_η` with `pV_z_to_ζ`
#     (and `pVI_η_to_z` with `pV_ζ_to_z`).  Verified bite (2 RED of
#     42): ET.1.5 — the η-frame problem's `zspan[1]` and `y0[2]` no
#     longer equal `pVI_z_to_η(z, u, up)`'s η and vp components (the
#     `to_frame` closure is now `pV_z_to_ζ`, returning `log z` not
#     `log log z`, and `z·up` not `z·log(z)·up`).
#
# Restoration: all three mutations restored before commit; GREEN at
# the new count.
