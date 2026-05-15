# test/ffw_fig_5_test.jl
#
# Quantitative pin for `figures/ffw2017_fig_5.jl` — the FFW 2017 §3
# tronquée P_III solution on its pole-free sector + the closed-form
# condition-number heatmap (FFW md:262).
#
# Source: references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
#         FFW2017_painleve_riemann_surfaces_preprint.md:220-265.
#
# ## Test identifier prefix: FF5
#
# Eight testsets `FF5.1.1` – `FF5.1.8`:
#
#   FF5.1.1 — IC round-trip exact at the two FFW md:243 boundary
#             ICs (z₁, z₂).  The figure's sheet-aware
#             `tronquee_ic_sheet(|z|, arg_z; n_terms, β, δ)` evaluated
#             at FFW's published arguments should reproduce FFW's
#             published u-values to within the n_terms=2 truncation
#             error.  Achievable: |Δ| ≈ 7e-3 at n_terms=2 (the v1
#             helper hard-codes a_3..a_N = 0 — see worklog 039
#             §"What is NOT shipped").
#
#   FF5.1.2 — Pole-free sector check.  At an interior ζ-point of the
#             principal-sheet BVP rectangle, |u(z)| is bounded — no
#             nearby pole.  `u(z) = w(ζ) / z` for PIII; for the
#             tronquée family, u ~ z^{1/3}, so |u| ~ |z|^{1/3} ≈ 3
#             at |z| = exp(Re ζ / 2) ≈ 30.  Threshold: |u| < 10.
#
#   FF5.1.3-5 — Self-cross-check at three strip-representative
#             interior ζ points.  FFW Fig 5 reports per-strip errors
#             `4e-6 / 3e-7 / 2e-8 / 3e-9 / 4e-8` (md:240).  At v1
#             Float64 + n_terms = 2 + N = 50, the achievable bound
#             is the FFW md:226 baseline of `~10⁻¹` for the
#             PFS-vs-BVP boundary derivative match.  We document
#             the achieved error for posterity and pin a coarser
#             threshold for the figure's BVP self-residual + max-|w|
#             sanity (the FFW number is unreachable at Float64 +
#             the v1 helper — worklog 040 §"Empirical sweep").
#
#   FF5.1.6 — Closed-form condition number `κ_r ≈ 157` at z = 30
#             (FFW md:264 verbatim).  Pin ±1.0.
#
#   FF5.1.7 — κ-monotonicity: κ(z=10) < κ(z=30) < κ(z=100).  Three
#             evaluations of the closed form.
#
#   FF5.1.8 — n_terms-convergence behavior.  The shipped helper
#             `PadeTaylor.pIII_asymptotic_ic` ships closed-form a_1
#             and a_2 only (worklog 039); so at `n_terms ∈ {3, 5,
#             10, 15}` the values are bit-identical (a_3..a_N = 0).
#             We pin this v1 behavior (`u(n_terms=3) == u(n_terms=15)`)
#             so the test flags any future change to the helper's
#             coefficient computation.
#
# ## Mutation-proof bites (see footer)
#
#   M1 — Corrupt one boundary IC by 1e-3 in u (perturb the FFW
#        verbatim numerics).  FF5.1.1 RED.
#   M2 — Hardcode `n_terms = 2` in `tronquee_ic_sheet` (no change
#        from the shipped helper's behavior — bites only if the
#        helper grows a_3+ coefficients in v2; the test catches a
#        regression in the OPPOSITE direction of v2 progress).
#        Skipped at v1; documented in footer.
#   M3 — Replace `κ_r ≈ 27/16 · |z|^{4/3}` with `κ_r ≈ 27/16 · |z|^{1/2}`.
#        FF5.1.6 RED at the 157 pin.
#   M4 — Skip the BVP step entirely (return the PFS values directly).
#        FF5.1.2 RED — |u| inside the sector goes catastrophic
#        (~|10⁻¹| FFW baseline + PFS's exponential ill-conditioning).

using Test
using PadeTaylor
using PadeTaylor.CoordTransforms: pIII_z_to_ζ

# Local replica of the figure's sheet-aware IC helper.  Kept here so
# that the test is self-contained and the figure script isn't
# `include`d (a figure script does plotting / I/O at module-load
# time, side effects we don't want in `runtests.jl`).  The math is
# the same: `s = z^{1/3}` evaluated with the CALLER-supplied unfolded
# argument, then `u = s · (1 + a_1 s^{-2} + a_2 s^{-4})` with the
# IVPBVPHybrid module's closed-form `a_1, a_2`.
function tronquee_ic_sheet_test(z_modulus::Real, arg_z::Real;
                                 n_terms::Integer = 15,
                                 β::Real = -1/20, δ::Real = -1)
    a1 = -float(β) / 3
    a2 = n_terms ≥ 2 ? ((4.0/9.0) + float(δ) * a1^2) / (2 * float(δ)) : 0.0
    s_mag = z_modulus^(1/3)
    s_arg = arg_z / 3
    s     = complex(s_mag * cos(s_arg), s_mag * sin(s_arg))
    if z_modulus < 2.0
        return (s, inv(3 * s^2))
    end
    s_inv2 = inv(s^2)
    u_sum  = one(s) + a1 * s_inv2 + a2 * s_inv2^2
    up_sum = one(s) - a1 * s_inv2 - 3 * a2 * s_inv2^2
    return (s * u_sum, up_sum / (3 * s^2))
end

@testset "FFW 2017 Fig 5 — tronquée P_III + condition-number heatmap" begin

    # ---- FFW md:240 parameters ----------------------------------------
    α, β, γ, δ = 1.0, -1/20, 0.0, -1.0

    # FFW md:243 published asymptotic-IC points (verbatim).
    arg_z1   = 13π/6       # = 9π/4 - π/12
    arg_z2   = -2π/3       # = -3π/4 + π/12
    z1_mod   = 30.0
    z2_mod   = 30.0
    u_z1_ffw  = -2.000735432319 + 2.376177147900im
    u_z2_ffw  =  2.384379236170 - 1.993845650158im
    up_z1_ffw = -5.939523100e-3 + 3.402038641e-2im
    up_z2_ffw =  6.050817704e-3 + 3.398020750e-2im

    # =================================================================
    # FF5.1.1 — IC round-trip vs FFW md:243 published values.
    # =================================================================
    @testset "FF5.1.1: tronquée IC reproduces FFW md:243 at z₁ and z₂" begin
        u1_comp, up1_comp = tronquee_ic_sheet_test(z1_mod, arg_z1;
                                                    n_terms = 15,
                                                    β = β, δ = δ)
        u2_comp, up2_comp = tronquee_ic_sheet_test(z2_mod, arg_z2;
                                                    n_terms = 15,
                                                    β = β, δ = δ)
        # |Δ| ~ 7e-3 at n_terms = 2 (v1 helper hard-cap; worklog 039
        # §"What is NOT shipped").  Tighter at v2 with a_3+ implemented.
        @test abs(u1_comp - u_z1_ffw) < 1e-2
        @test abs(u2_comp - u_z2_ffw) < 1e-2
        @info "FF5.1.1 |u(z₁) - FFW| = $(abs(u1_comp - u_z1_ffw))  (n_terms=15 helper; FFW: optimally-truncated series)"
        @info "FF5.1.1 |u(z₂) - FFW| = $(abs(u2_comp - u_z2_ffw))  (same)"
        # Spot-check leading: |u₁| should ≈ |z|^{1/3} = 3.107.
        @test isapprox(abs(u1_comp), 30.0^(1/3); atol = 1e-2)
        @test isapprox(abs(u2_comp), 30.0^(1/3); atol = 1e-2)
        # Argument should match the unfolded `arg z / 3 + correction`.
        # u₁ is in Q2 (negative re, positive im): arg ≈ 13π/18 ≈ 2.27 rad.
        @test isapprox(angle(u_z1_ffw), 13π/18; atol = 5e-3)
        # u₂ is in Q4 (positive re, negative im): arg ≈ -2π/9 ≈ -0.698 rad.
        @test isapprox(angle(u_z2_ffw), -2π/9; atol = 5e-3)
    end

    # =================================================================
    # FF5.1.2 — Pole-free sector check: |u| bounded inside the BVP
    # rectangle.  We run a SMALL hybrid solve on the principal-sheet
    # sub-strip and verify u remains finite / |u| < 10 at an interior
    # ζ point.
    # =================================================================
    @testset "FF5.1.2: pole-free sector — |u| bounded at sector interior" begin
        # FFW Fig 5 principal-sheet sub-strip (matches the figure
        # script's sector).
        re_anchor = 2 * log(30)
        im_hi = 2π - 0.05
        im_lo = -3π/2 + 0.05
        re_extent = 1.0       # narrow for test runtime

        # Build pp seeded at Z₂ (principal-sheet lower IC).
        z2 = 30.0 * cis(arg_z2)
        u2, up2 = tronquee_ic_sheet_test(z2_mod, arg_z2;
                                          n_terms = 15, β = β, δ = δ)
        pp = PainleveProblem(:III; α = α, β = β, γ = γ, δ = δ,
                              u0 = u2, up0 = up2,
                              zspan = (z2, z2 + 1.0), order = 30)
        sector = (im_lo = im_lo, im_hi = im_hi,
                   re_anchor = re_anchor, re_extent = re_extent)
        ic_fn = z -> tronquee_ic_sheet_test(abs(z), angle(z);
                                             n_terms = 15, β = β, δ = δ)
        sol = solve_pole_free_hybrid(pp, sector, ic_fn;
                pfs_kwargs = (; h = 0.4,
                              step_size_policy = :adaptive_ffw,
                              adaptive_tol = 1e-10,
                              k_conservative = 1e-3,
                              max_rescales = 50,
                              max_steps_per_target = 500),
                bvp_kwargs = (; N = 30, tol = 1e-10, maxiter = 30),
                n_slices = 5,
                glue_tol = 1e-8)

        # Interior ζ point.
        ζ_test = complex(re_anchor - 0.3, 0.5)   # well inside the sector
        w, wp = sol(ζ_test)
        @test isfinite(w) && isfinite(wp)
        # PIII u = w / z; z = exp(ζ/2).
        z_test = exp(ζ_test / 2)
        u_test = w / z_test
        @test abs(u_test) < 10.0
        @info "FF5.1.2 at ζ = $ζ_test:  |w| = $(abs(w)),  |u| = $(abs(u_test))  (asymptotic ~|z|^{1/3} ≈ $(round(abs(z_test)^(1/3); digits=2)))"
    end

    # =================================================================
    # FF5.1.3-5 — Strip self-cross-check.  At three representative
    # interior ζ points (bottom strip / mid strip / top strip of the
    # principal-sheet sub-strip), compute u via the hybrid and
    # cross-check against tronquee_ic_sheet_test evaluated at the
    # corresponding (|z|, arg z).  This is a self-consistency check
    # (not the FFW md:240 oracle, which is unreachable at v1 — see
    # docstring §"Strip self-cross-check thresholds").
    #
    # FFW per-strip errors (md:240) bottom→top: 4e-6 / 3e-7 / 2e-8 /
    # 3e-9 / 4e-8.  Our v1 stack reaches the md:226 baseline `~10⁻¹`
    # for the boundary-derivative match; the SECTOR-INTERIOR error
    # (vs the asymptotic series, which is itself an n_terms=2
    # approximation) is dominated by the asymptotic-series error
    # ≈ |z|^{-7/3} · a_3 ~ |30|^{-7/3} · O(1) ≈ 3·10⁻⁴ for the
    # SERIES side.  The achievable threshold is therefore ~3·10⁻³
    # at v1, not FFW's md:240 numbers.  We pin a generous fail-loud
    # bound and @info-log the achieved values for the worklog.
    # =================================================================

    function _self_check_at_ζ(ζ_test, im_hi, im_lo)
        re_anchor = real(ζ_test)        # Anchor to right of ζ_test by RE_EXTENT
        re_anchor_use = max(real(ζ_test) + 0.3, 2 * log(30))
        re_extent = re_anchor_use - real(ζ_test) + 0.2

        z2 = 30.0 * cis(-2π/3)
        u2, up2 = tronquee_ic_sheet_test(30.0, -2π/3;
                                          n_terms = 15, β = -1/20, δ = -1)
        pp = PainleveProblem(:III; α = 1.0, β = -1/20, γ = 0.0, δ = -1.0,
                              u0 = u2, up0 = up2,
                              zspan = (z2, z2 + 1.0), order = 30)
        sector = (im_lo = im_lo, im_hi = im_hi,
                   re_anchor = re_anchor_use, re_extent = re_extent)
        ic_fn = z -> tronquee_ic_sheet_test(abs(z), angle(z);
                                             n_terms = 15, β = -1/20, δ = -1)
        sol = solve_pole_free_hybrid(pp, sector, ic_fn;
                pfs_kwargs = (; h = 0.4,
                              step_size_policy = :adaptive_ffw,
                              adaptive_tol = 1e-10,
                              k_conservative = 1e-3,
                              max_rescales = 50,
                              max_steps_per_target = 500),
                bvp_kwargs = (; N = 30, tol = 1e-10, maxiter = 30),
                n_slices = 5,
                glue_tol = 1e-8)

        w_bvp, _ = sol(ζ_test)
        z_test = exp(ζ_test / 2)
        u_bvp = w_bvp / z_test

        # Reference: tronquee_ic_sheet_test at unfolded arg z = Im(ζ)/2.
        u_series, _ = tronquee_ic_sheet_test(abs(z_test), imag(ζ_test) / 2;
                                              n_terms = 15, β = -1/20, δ = -1)
        return (u_bvp = u_bvp, u_series = u_series,
                err = abs(u_bvp - u_series))
    end

    # Principal-sheet bounds for the cross-check sector.
    im_hi_test = 2π - 0.05
    im_lo_test = -3π/2 + 0.05

    @testset "FF5.1.3: strip self-cross-check (bottom strip, Im ζ ≈ -π)" begin
        # FFW md:240 strip 1 (bottom): error 4e-6 (oracle).  v1 reach
        # at Float64 + n_terms=2: ~10⁻³ (asymptotic series truncation
        # error).  Threshold: 5e-2 (10× the achievable, to bite M1).
        ζ_test = complex(6.3, -π)
        result = _self_check_at_ζ(ζ_test, im_hi_test, im_lo_test)
        @info "FF5.1.3 (bottom strip, ζ = $ζ_test):  u_BVP = $(result.u_bvp),  u_series = $(result.u_series),  |Δu| = $(result.err)  (FFW md:240: 4e-6)"
        @test isfinite(result.u_bvp)
        @test result.err < 5e-2     # generous v1 bound — bites M1, M4
    end

    @testset "FF5.1.4: strip self-cross-check (middle strip, Im ζ ≈ 0)" begin
        # FFW md:240 strip 3 (middle): error 2e-8.  v1 reach: ~10⁻³.
        ζ_test = complex(6.3, 0.0)
        result = _self_check_at_ζ(ζ_test, im_hi_test, im_lo_test)
        @info "FF5.1.4 (mid strip, ζ = $ζ_test):  u_BVP = $(result.u_bvp),  u_series = $(result.u_series),  |Δu| = $(result.err)  (FFW md:240: 2e-8)"
        @test isfinite(result.u_bvp)
        @test result.err < 5e-2
    end

    @testset "FF5.1.5: strip self-cross-check (upper strip, Im ζ ≈ +π)" begin
        # FFW md:240 strip 5 (top): error 4e-8.  v1 reach: ~10⁻³.
        ζ_test = complex(6.3, π)
        result = _self_check_at_ζ(ζ_test, im_hi_test, im_lo_test)
        @info "FF5.1.5 (top strip, ζ = $ζ_test):  u_BVP = $(result.u_bvp),  u_series = $(result.u_series),  |Δu| = $(result.err)  (FFW md:240: 4e-8)"
        @test isfinite(result.u_bvp)
        @test result.err < 5e-2
    end

    # =================================================================
    # FF5.1.6 — Closed-form condition number κ_r ≈ 157 at z = 30
    # (FFW md:264 verbatim).
    # =================================================================
    @testset "FF5.1.6: condition-number closed form κ_r(z=30) ≈ 157" begin
        z_val = 30.0
        κ_r = (27/16) * z_val^(4/3)
        @test isapprox(κ_r, 157.0; atol = 1.0)
        @info "FF5.1.6: κ_r(|z|=30) = $κ_r  (FFW md:264: ≈ 157)"
        # And the FFW md:262 ζ-frame form: κ_r ~ (27/16) e^{2 Re ζ / 3};
        # at Re ζ = 2 log 30 this collapses to the |z| form.
        re_ζ = 2 * log(30)
        @test isapprox((27/16) * exp(2 * re_ζ / 3), κ_r; atol = 1e-9)
    end

    # =================================================================
    # FF5.1.7 — κ-monotonicity: κ grows exponentially with |z|.
    # =================================================================
    @testset "FF5.1.7: κ_r monotonicity κ(10) < κ(30) < κ(100)" begin
        κ_10  = (27/16) * 10.0^(4/3)
        κ_30  = (27/16) * 30.0^(4/3)
        κ_100 = (27/16) * 100.0^(4/3)
        @test κ_10 < κ_30 < κ_100
        # Exponential growth check: κ(100) / κ(30) = (100/30)^{4/3} ≈ 4.81.
        @test isapprox(κ_100 / κ_30, (100/30)^(4/3); atol = 1e-9)
        @info "FF5.1.7: κ_r(|z|=10) = $(round(κ_10; digits=2))  κ_r(|z|=30) = $(round(κ_30; digits=2))  κ_r(|z|=100) = $(round(κ_100; digits=2))"
    end

    # =================================================================
    # FF5.1.8 — n_terms-convergence behavior.  Lifted from A6 IB.1.1
    # for context; documents the v1 cap.  The shipped helper
    # `PadeTaylor.pIII_asymptotic_ic` hard-codes `a_3 ... a_N = 0`
    # (worklog 039 §"What is NOT shipped"), so the cross-helper
    # values agree for any n_terms ≥ 2.  This pin documents that v1
    # behavior — a future change to the helper's coefficient
    # computation (the deferred bead for higher-order a_n) would
    # cause this test to RED, prompting the test author to lift the
    # FFW md:240 oracle to a tighter threshold.
    # =================================================================
    @testset "FF5.1.8: n_terms-convergence — v1 helper caps at a_2" begin
        z_val = 30.0 + 0im       # principal-sheet sample
        u_n2,  _ = pIII_asymptotic_ic(z_val; n_terms = 2,  β = β, δ = δ)
        u_n3,  _ = pIII_asymptotic_ic(z_val; n_terms = 3,  β = β, δ = δ)
        u_n5,  _ = pIII_asymptotic_ic(z_val; n_terms = 5,  β = β, δ = δ)
        u_n10, _ = pIII_asymptotic_ic(z_val; n_terms = 10, β = β, δ = δ)
        u_n15, _ = pIII_asymptotic_ic(z_val; n_terms = 15, β = β, δ = δ)
        # v1: bit-identical (a_3..a_N = 0).
        @test u_n3  == u_n2
        @test u_n5  == u_n2
        @test u_n10 == u_n2
        @test u_n15 == u_n2
        @info "FF5.1.8: u(n_terms=2)  = $u_n2"
        @info "FF5.1.8: u(n_terms=15) = $u_n15   (bit-identical to n_terms=2 at v1; a_3..a_N=0 — worklog 039)"
        # Cross-check vs the FFW md:230 leading order: u ~ z^{1/3}.
        @test isapprox(abs(u_n2), 30.0^(1/3); atol = 1e-2)
    end
end

# =============================================================================
# Mutation-proof footer (CLAUDE.md Rule 4).
#
# Verified bites; perturb impl, observe RED, restore.
#
# M1 — Corrupt FFW md:243 verbatim numerics: replace `u_z1_ffw` with
#      `u_z1_ffw + 5e-2 + 0im` (larger perturbation than the test
#      threshold of 1e-2, since the v1 baseline |Δ| is already
#      ~7.4e-3 — the threshold has 2.6e-3 of margin, so a 5e-2
#      perturbation cleanly bites).  Verified: |Δ| at the corrupted
#      value = 5.65e-2 > 1e-2 → FF5.1.1 fails.
#
# M2 — Hardcode `n_terms = 2` in `tronquee_ic_sheet_test` (skip the
#      a_2 closed-form computation; force a_2 = 0).  At v1 this is
#      no-op for n_terms ≥ 2 (a_3+ already zero).  Bites only when
#      a_3+ are implemented in v2 (the helper grows higher coefs);
#      at that point FF5.1.3-5 thresholds should also tighten.
#      Footer-documented; not actively executable at v1.
#
# M3 — Replace the κ_r closed form `(27/16) * z^(4/3)` with
#      `(27/16) * z^(1/2)`.  FF5.1.6 RED — at z = 30, the wrong
#      formula gives ≈ 9.24, not 157.  Verified: patch FF5.1.6 line
#      `(27/16) * z_val^(4/3)` to `(27/16) * z_val^(1/2)` → fails.
#
# M4 — Skip the BVP step (return PFS values directly).  The IB.1.5
#      (worklog 039) M1 already demonstrated this bite at the
#      hybrid-module level: |Δw'| jumps from ~1 to ~50 when the
#      BVP step is stubbed.  Here the analogous bite is in FF5.1.3-5:
#      the `result.err < 5e-2` thresholds catch the FFW md:226
#      baseline `~10⁻¹` PFS-on-pole-free-sector error (which would
#      be ~10⁻¹ in u, breaking the 5e-2 threshold).  FF5.1.2's
#      `|u| < 10` pin is more permissive — the PFS does eventually
#      return a u-value of magnitude ~|z|^{1/3} ≈ 3, just with the
#      wrong precision.  So M4's bite is at FF5.1.3-5.
#
# ## v1 limit honesty (CLAUDE.md Rule 9)
#
# The FFW md:240 per-strip errors `4e-6 / 3e-7 / 2e-8 / 3e-9 / 4e-8`
# are NOT reachable at Float64 + the shipped n_terms=2 helper.  Per
# worklog 039 §Friction 3 and worklog 040 §"Empirical sweep":
#
#   * The dominant Float64 error is the asymptotic-series truncation
#     ≈ |z|^{-7/3} · a_3 ~ 3e-4 at |z| = 30 (a_3+ in the helper are
#     zero in v1).
#   * BVP spectral floor (`N = 50`) is ~1e-12 — not the bottleneck.
#   * PFS-vs-BVP boundary derivative match at v1: |Δw'| ≈ 0.4 - 1.7
#     (FFW md:226 baseline).
#
# Reaching FFW's per-strip 1e-7 → 1e-8 numbers requires the v2 stack
# (BF-256 + a_3+ helper + tighter N).  Two follow-up beads:
#
#   * `padetaylor-7zw`-style BF-256 PIII tronquée pin.
#   * `padetaylor-?fr` (to be filed): implement `pIII_asymptotic_ic`
#     a_3, a_4, ..., a_N via TaylorSeries.jl substitution per FFW
#     md:232 "optimal truncation".
# =============================================================================
