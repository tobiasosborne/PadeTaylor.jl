# test/problems_test.jl -- Phase 6 tests.
#
# The Padé-vs-Taylor pole-bridge demonstration: ONE Padé built at z = 0
# with h_max = 1.5 spans the lattice pole of u(z) = ℘(z + c₁; 0, c₂)
# at z = 1.  In the rescaled variable t = z / h_max the pole sits at
# t = 1/1.5 ≈ 0.667, strictly INSIDE the segment [0, 1].  The same Padé,
# evaluated at t = z/1.5 for z ∈ [0, 1.5], must therefore give correct
# values both BEFORE the pole (z = 0.5, 0.95) and AFTER (z = 1.05, 1.4)
# — bridging the pole is the central claim of the analytic-continuation
# advantage of Padé over plain Taylor.
#
# Reference: docs/worklog/004-phase-6-pivot.md (full spec, oracle plan,
# v1/v2 acceptance scope).  Oracles cross-checked in
# external/probes/problems-oracle/verify.jl.
#
# v1 acceptance: tests 6.1.1–6.1.6 GREEN at the listed tolerances.  v2
# (FW Table 5.1 long-range integration to z = 30) is deferred — see
# bead padetaylor-8cr.

using Test
using PadeTaylor

include(joinpath(@__DIR__, "_oracle_problems.jl"))

@testset "Problems (Phase 6): Padé bridges the lattice pole at z=1" begin
    fW(z, u, up) = 6 * u^2

    @testset "6.1.1  sol(0.5) ~ WeierstrassP(0.5) to 1e-13 [3-source]" begin
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
        sol = solve_pade(prob; h = 1.5)
        u, up = sol(0.5)
        @test isapprox(u,  u_at_0_5;  rtol = 1e-13)
        @test isapprox(up, up_at_0_5; rtol = 1e-13)
    end

    @testset "6.1.2  sol(0.95) ~ WeierstrassP(0.95) to 1e-9 [near-pole, 3-source]" begin
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
        sol = solve_pade(prob; h = 1.5)
        u, up = sol(0.95)
        @test isapprox(u,  u_at_0_95;  rtol = 1e-9)
        @test isapprox(up, up_at_0_95; rtol = 1e-9)
    end

    @testset "6.1.3  sol(1.05) ~ WeierstrassP(1.05) to 1e-9 -- Padé bridges pole" begin
        # Headline: same Padé, evaluated at t = 1.05/1.5 = 0.7, i.e. just
        # PAST the pole at t = 0.667.  Plain Taylor diverges here (test 6.1.5);
        # Padé tracks the closed-form to 1e-9 because the denominator zero
        # at t ≈ 0.667 captures the lattice pole structure.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
        sol = solve_pade(prob; h = 1.5)
        u, up = sol(1.05)
        @test isapprox(u,  u_at_1_05;  rtol = 1e-9)
        @test isapprox(up, up_at_1_05; rtol = 1e-9)
    end

    @testset "6.1.4  sol(1.4) ~ WeierstrassP(1.4) to 1e-7 -- further past pole" begin
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
        sol = solve_pade(prob; h = 1.5)
        u, up = sol(1.4)
        @test isapprox(u,  u_at_1_4;  rtol = 1e-7)
        @test isapprox(up, up_at_1_4; rtol = 1e-7)
    end

    @testset "6.1.5  Padé wins, plain Taylor diverges -- the headline demo" begin
        # Same Taylor coefficients on both sides; the only difference is
        # whether we Padé-convert or truncate.  Past the natural radius of
        # convergence at z = 1, plain Taylor explodes by orders of magnitude;
        # Padé continues to track the closed-form ℘.
        prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
        sol = solve_pade(prob; h = 1.5)
        u_pade, _ = sol(1.05)

        coefs_u = taylor_coefficients_2nd(fW, 0.0, u_0_FW, up_0_FW, 30)
        u_taylor = PadeTaylor.taylor_eval(coefs_u, 1.05)

        # Padé wins: rel err < 1e-9.
        @test isapprox(u_pade, u_at_1_05; rtol = 1e-9)
        # Taylor loses: relative error > 0.1 (truncation diverges past
        # natural radius of convergence at z = 1).
        @test abs(u_taylor - u_at_1_05) / abs(u_at_1_05) > 0.1
    end

    @testset "6.1.6  BigFloat-256: sol(1.05) at high precision (order = 40)" begin
        # Why order = 40, not 30: empirical sweep at z = 1.05 (just past the
        # pole) shows the (15, 15) diagonal Padé from order = 30 has an
        # inherent rational-approximation-error floor of ~3.5e-10 at
        # t = 0.7 — at THAT order BF-256 adds no value over Float64 because
        # both saturate at the same Padé-error floor.  Bumping to order = 40
        # changes the picture: Float64 *degrades* to ~1.6e-8 (SVD ill-
        # conditioning of the 41-by-41 Toeplitz block), while BF-256
        # cleanly converges to ~5e-15 (the IC-precision floor for the
        # 16-digit FW ICs after amplification through the pole).  This
        # makes 6.1.6 a meaningful test of the BF-256 path: at this
        # order, BF-256 is REQUIRED, not optional.
        setprecision(BigFloat, 256) do
            u0_bf  = BigFloat("1.071822516416917")
            up0_bf = BigFloat("1.710337353176786")
            zstart = BigFloat("0")
            zend   = BigFloat(3) / BigFloat(2)   # exact 1.5

            prob = PadeTaylorProblem(fW, (u0_bf, up0_bf), (zstart, zend); order = 40)
            sol = solve_pade(prob; h = BigFloat(3) / BigFloat(2))

            u_ref  = parse(BigFloat, u_at_1_05_80dps_str)
            up_ref = parse(BigFloat, up_at_1_05_80dps_str)

            u_eval, up_eval = sol(BigFloat(105) / BigFloat(100))
            @test isapprox(u_eval,  u_ref;  rtol = BigFloat("1e-13"))
            @test isapprox(up_eval, up_ref; rtol = BigFloat("1e-13"))
        end
    end
end

# ===========================================================================
# 6.2  DESCENDING-span support (bug padetaylor-xhjw).
# Pre-fix, solve_pade with z_end < z_start returned a DEGENERATE one-node
# trajectory (the loop `while state.z < z_end` never entered) — a silent
# Rule-1 violation.  Now the driver steps with sign(z_end - z_start)·h_max and
# the callable is direction-agnostic.  Oracle: closed-form ℘ (no code shared
# with src/), test/_oracle_problems.jl §C, wolframscript-regenerated; verified
# ℘(-1;{0,2}) == u_0_FW.  Tolerances are the same spectral floors the ascending
# 6.1.x pins hit (not fitted): rtol 1e-11 (u) / 1e-9 (the smaller-magnitude u').
# ===========================================================================
@testset "Problems 6.2: descending-span support (xhjw)" begin
    fW(z, u, up) = 6 * u^2
    y0 = (u_0_FW, up_0_FW)

    # GM-1/EC-1/EC-2: single descending segment (0,-0.5), h_max=0.5.
    @testset "6.2.1  descending single segment: value + exact endpoint" begin
        prob = PadeTaylorProblem(fW, y0, (0.0, -0.5); order = 30)
        sol  = solve_pade(prob; h = 0.5)
        @test length(sol.h) == 1               # EC-1: one step, NOT degenerate
        @test sol.z == [0.0, -0.5]             # EC-2: no spurious residue step
        @test sol.z[end] == -0.5               # snap: endpoint EXACT
        u3, up3 = sol(-0.3)                     # GM-1: interior value vs ℘
        @test isapprox(u3,  u_at_m0_3;  rtol = 1e-11)
        @test isapprox(up3, up_at_m0_3; rtol = 1e-9)
        u5, up5 = sol(-0.5)                     # GM-2: descending endpoint vs ℘
        @test isapprox(u5,  u_at_m0_5;  rtol = 1e-11)
        @test isapprox(up5, up_at_m0_5; rtol = 1e-9)
        @test isapprox(sol(0.0)[1], u_0_FW; rtol = 1e-13)   # EC-3: z_start = IC
    end

    # GM-3: MULTI-segment descending — clamp (segment count) AND the callable's
    # direction-aware segment scan.  Two spans: a smooth one (value + clamp) and
    # a POLE-BRIDGE one where the scan is load-bearing — a point PAST the pole is
    # valued correctly ONLY by the pole-bridging segment; an ascending-only scan
    # would pick the wrong segment (off by ~0.06 at z=-1.8).
    @testset "6.2.2  descending multi-segment: clamp + scan (incl. pole-bridge)" begin
        # (a) smooth span (0,-1.0), h=0.2 ⇒ 5 segments (clamp: NOT one big step).
        sol = solve_pade(PadeTaylorProblem(fW, y0, (0.0, -1.0); order = 30); h = 0.2)
        @test length(sol.h) == 5
        @test sol.z[1] == 0.0 && sol.z[end] == -1.0
        @test issorted(sol.z; rev = true)            # strictly DECREASING
        @test isapprox(sol(-0.7)[1], u_at_m0_7; rtol = 1e-10)
        @test isapprox(sol(-0.3)[1], u_at_m0_3; rtol = 1e-10)
        # (b) POLE-BRIDGE span (0,-1.95), h=1.0 ⇒ 2 segments; the lattice pole at
        # z=-1.726 sits in segment 2.  sol(-1.8) is PAST the pole — only
        # segment 2's pole-bridging Padé values it right, so the scan MUST select
        # it (the assertion that makes the direction-aware scan load-bearing).
        solp = solve_pade(PadeTaylorProblem(fW, y0, (0.0, -1.95); order = 30); h = 1.0)
        @test length(solp.h) == 2
        @test solp.z[end] == -1.95
        @test isapprox(solp(-1.8)[1], u_at_m1_8; rtol = 1e-7)   # scan-sensitive
        @test isapprox(solp(-0.5)[1], u_at_m0_5; rtol = 1e-10)  # segment 1, smooth
    end

    # EC-4: just-outside on BOTH sides throws the reframed [lo,hi] window guard.
    @testset "6.2.3  window guard rejects just-outside on both ends" begin
        prob = PadeTaylorProblem(fW, y0, (0.0, -0.5); order = 30)
        sol  = solve_pade(prob; h = 0.5)
        @test_throws DomainError sol(1e-9)         # above hi = 0.0
        @test_throws DomainError sol(-0.5 - 1e-9)  # below lo = -0.5
    end

    # GM-4: forward∘reverse round-trip recovers the IC (no external oracle) —
    # the krgy.3 MM.5 reversibility invariant at the solve_pade level.
    @testset "6.2.4  forward∘reverse round-trip recovers the IC" begin
        solA = solve_pade(PadeTaylorProblem(fW, y0, (0.0, 0.5); order = 30); h = 1.5)
        uE, upE = solA(0.5)
        solB = solve_pade(PadeTaylorProblem(fW, (uE, upE), (0.5, 0.0); order = 30); h = 1.5)
        u0, up0 = solB(0.0)
        @test isapprox(u0,  u_0_FW;  rtol = 1e-11)
        @test isapprox(up0, up_0_FW; rtol = 1e-11)
    end

    # GM-5: BigFloat-256 descending round-trip — proves the fix is precision-
    # clean (sign(), the dir·min clamp, and the state.z = z_end snap are all
    # generic-T; no Float64 fallback in the descending path).
    @testset "6.2.5  BigFloat-256 descending round-trip is precision-clean" begin
        setprecision(BigFloat, 256) do
            y0b = (BigFloat("1.071822516416917"), BigFloat("1.710337353176786"))
            h   = BigFloat(3) / BigFloat(2)
            half = BigFloat(1) / 2                       # exact 0.5 in BigFloat
            solA = solve_pade(PadeTaylorProblem(fW, y0b, (BigFloat(0), half);
                                                order = 40); h = h)
            @test solA.z[end] == half                    # snap exact in BigFloat
            uE, upE = solA(half)
            solB = solve_pade(PadeTaylorProblem(fW, (uE, upE), (half, BigFloat(0));
                                                order = 40); h = h)
            u0, up0 = solB(BigFloat(0))
            @test isapprox(u0,  y0b[1]; rtol = BigFloat("1e-13"))
            @test isapprox(up0, y0b[2]; rtol = BigFloat("1e-13"))
        end
    end

    # EC-6: a zero-length span stays rejected at CONSTRUCTION — the descending
    # support opened no new empty-span path.
    @testset "6.2.6  z_end == z_start still rejected at construction" begin
        @test_throws ArgumentError PadeTaylorProblem(fW, y0, (0.3, 0.3); order = 30)
    end
end

# 6.1.7  Mutation-proof procedure (verified manually before commit; see
# worklog 004 §"Mutation-proof procedure for the next agent"):
#
#   Mutation A  --  in `pade_step_with_pade!` (or `solve_pade`), replace
#     the `robust_pade(c̃, m, n)` call with one that returns the Taylor-
#     truncation Padé (numerator = c̃, denominator = [1.0]).
#     Expected: 6.1.3 + 6.1.4 RED (Taylor diverges past the pole at t=0.667).
#
#   Mutation B  --  in `(sol::PadeTaylorSolution)(z)`, return the
#     segment-start state `sol.y[k]` instead of evaluating the Padé at
#     the rescaled t.
#     Expected: 6.1.1 RED  (z = 0.5 returns u(0) ≈ 1.072 not u(0.5) ≈ 4.004).
#
# 6.2.7  Mutation-proof of the DESCENDING-span fix (bug padetaylor-xhjw,
#   2026-06-13) — EXECUTED on src/Problems.jl, each restored byte-clean
#   (`git diff src/Problems.jl` then showed only the fix).  Four independent
#   components, each with a distinct RED:
#     MT-1 LOOP   — revert `while dir*(z_end-state.z) > zero(T)` to
#       `while state.z < z_end`: the descending loop never enters ⇒ degenerate
#       one-node trajectory ⇒ 6.2.1 (steps==1) + value evals RED; 6.1.x GREEN.
#     MT-2 CLAMP  — revert `dir*min(h_max_T, abs(gap))` to `min(h_max_T, gap)`:
#       for a descending span `min` picks the (negative) full gap ⇒ ONE giant
#       step ⇒ 6.2.2 `length(sol.h)==5` / `==2` RED.
#     MT-3 SCAN   — revert the callable scan `dir*(z_T - sol.z[k+1]) > zero(T)`
#       to `z_T > sol.z[k+1]`: descending multi-segment picks the WRONG segment.
#       This is INVISIBLE on the smooth pole-free trajectory (any segment's
#       (15,15) Padé extrapolates ℘ to ~2e-13 even at t≈9), so the load-bearing
#       assertion is the POLE-BRIDGE 6.2.2(b): sol(-1.8) past the z=-1.726 pole
#       is valued ONLY by segment 2's pole-bridging Padé (segment 1 is off by
#       ~0.06) ⇒ RED.
#     MT-4 GUARD  — revert `lo,hi = minmax(sol.z[1], sol.z[end])` to
#       `lo,hi = sol.z[1], sol.z[end]`: for descending sol.z, lo>hi ⇒ the live
#       bug (sol(-0.3) throws "outside window") ⇒ 6.2.1/6.2.3 RED.
#   NOTE (no snap test): an earlier draft snapped `state.z = z_end` on the final
#   step.  A sweep of 570 (z_end, h_max) spans found ZERO where it changed
#   behaviour — the clamp's `h_step == gap` already lands exactly
#   (`x + (z_end - x) == z_end`) — so the snap was REMOVED as verified-redundant
#   (Rule 9: no untestable code).  Endpoint exactness is pinned by the
#   `sol.z[end] == z_end` assertions, which the MT-2 clamp mutation reddens.
#   SIBLING: src/VectorProblems.jl:262 has the identical forward-only bug,
#   tracked as padetaylor-x0p0 (its adaptive Jorba-Zou stepper needs a distinct
#   fix + battery).
