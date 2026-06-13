"""
V3b tests for `PadeTaylor.VectorProblems` (bead `padetaylor-0ln.9`;
v0.2 plan row V3b).

`vector_solve_pade(prob; h)` is the top-level driver for first-order
vector ODEs `y' = f(z, y)`, `y ∈ ℂ^d`.  It loops `vector_pade_step!`
over a fixed step `h`, stores a per-segment shared-`Q` approximant, and
exposes a dense-output callable `sol(z)`.

These tests (`VP.*`) are written RED-first per CLAUDE.md Rule 4; each
asserts an invariant against a known-correct value (Rule 5):

    VP.1.1  d=1 reduction        — a 1-component solve reproduces the
                                   explicit scalar Padé–Taylor pipeline
    VP.1.2  closed-form system   — harmonic system → [cos z, −sin z]
                                   at every breakpoint AND via sol(z)
    VP.1.3  meromorphic system   — d=2 system with a genuine shared
                                   pole; integrate a smooth stretch
    VP.1.4  dense callable       — sol(z) interior + endpoints + range
    VP.1.5  failure modes        — h ≤ 0, order < 2, coincident zspan,
                                   max_steps exceeded
    VP.1.6  mutation-proof       — recorded at end of file
    VP.2.*  descending span      — bug `padetaylor-x0p0`: the dir-aware driver
                                   + callable integrate a DESCENDING span
                                   (`z_end < z_start`) instead of silently
                                   returning a degenerate one-node trajectory;
                                   oracle = the equianharmonic ℘ companion
    VP.2.7  mutation-proof        — the x0p0 fix (MT-1..MT-4), recorded at end

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/vector_problems_test.jl`) and under
`runtests.jl`.

NOTE — V3a found that for *entire* (pole-free) systems the shared `Q`
carries a small fixed fit residual (~4e-12 Float64) that does not shrink
with order or precision.  The VP.1.2 tolerances carry that headroom; it
is documented, not a bug.
"""

using Test
using PadeTaylor
using PadeTaylor.VectorProblems:
    VectorPadeTaylorProblem, VectorPadeTaylorSolution, vector_solve_pade
using PadeTaylor.Coefficients: taylor_coefficients_1st
using PadeTaylor.RobustPade:   robust_pade

# -----------------------------------------------------------------------------
# Test-local helpers (inline oracles — kept independent of the impl).
# -----------------------------------------------------------------------------

# Evaluate the polynomial `c` (low-to-high coefficients) at `t` by Horner.
_horner(c, t) = begin
    s = zero(eltype(c)) * t
    for k in length(c):-1:1
        s = s * t + c[k]
    end
    s
end

# Rescale c̃_k = h^k · c_k (the variable substitution h' = h·t).
_rescale(c, h) = [c[k+1] * h^k for k in 0:(length(c)-1)]

# The explicit scalar Padé pipeline, built inline for the VP.1.1 oracle:
# one fixed-h step — Taylor jet → h^k rescale → robust_pade(:svd) →
# evaluate at t = 1.
function _scalar_pade_step(f, z0, y0, order, h)
    jet = taylor_coefficients_1st(f, z0, y0, order)
    scaled = _rescale(jet, h)
    m = order ÷ 2
    P = robust_pade(scaled, m, m; method = :svd)
    return _horner(P.a, one(h)) / _horner(P.b, one(h))
end

@testset "VectorProblems (V3b)" begin

    @testset "VP.1.1 d=1 reduces to the explicit scalar pipeline" begin
        # A 1-component problem solved by `vector_solve_pade` must
        # reproduce the explicit scalar Padé–Taylor pipeline step by
        # step.  Scalar ODE: y' = z² + y², y(0) = 0 (FW 2011 §2.2.1),
        # order 30.  Integrate [0, 0.6] in two h = 0.3 steps.
        scalar_f = (z, y) -> z^2 + y^2
        vec_f    = (z, y) -> [z^2 + y[1]^2]
        order, h = 30, 0.3

        prob = VectorPadeTaylorProblem(vec_f, [0.0], (0.0, 0.6); order)
        sol  = vector_solve_pade(prob; h)

        # Replay the scalar pipeline segment by segment.
        y_ref = 0.0
        z_ref = 0.0
        @test sol.z[1] ≈ 0.0
        @test sol.y[1][1] ≈ 0.0
        for k in 1:length(sol.h)
            y_ref = _scalar_pade_step(scalar_f, z_ref, y_ref, order, h)
            z_ref += h
            @test sol.z[k+1]   ≈ z_ref
            @test sol.y[k+1][1] ≈ y_ref atol = 1e-12
        end
        @test sol.z[end] ≈ 0.6
    end

    @testset "VP.1.2 closed-form harmonic system over an interval" begin
        # y₁' = y₂, y₂' = −y₁, y0 = [1, 0] ⇒ y(z) = [cos z, −sin z].
        # Both components entire (no poles) — under the correct GGT (m,m) window
        # the shared Q reduces to its honest supported degree and carries a
        # method-set residual ~5e-9 (F64) / ~1.6e-17 (BF) for this pole-free
        # system (ADR-0027 / bead padetaylor-unk; cf. V3a VS.1.2).  Dispatch
        # (ADR-0028) will recover the tighter accuracy via the (m−1,m) cell.
        harm = (z, y) -> [y[2], -y[1]]
        order = 24

        prob = VectorPadeTaylorProblem(harm, [1.0, 0.0], (0.0, 2.0); order)
        sol  = vector_solve_pade(prob; h = 0.25)

        # Trajectory at every breakpoint.
        for k in 1:length(sol.z)
            z = sol.z[k]
            @test sol.y[k][1] ≈ cos(z)  atol = 1e-8
            @test sol.y[k][2] ≈ -sin(z) atol = 1e-8
        end
        @test sol.z[end] ≈ 2.0

        # Dense callable at several interior points.
        for z in (0.1, 0.37, 0.99, 1.5, 1.875)
            yz = sol(z)
            @test yz[1] ≈ cos(z)  atol = 1e-8
            @test yz[2] ≈ -sin(z) atol = 1e-8
        end

        # BigFloat (256-bit): far below Float64 reach.  The sin-jet
        # shared-Q residual is ~3.4e-21 (V3a VS.1.2); tolerance reflects.
        setprecision(BigFloat, 256) do
            probb = VectorPadeTaylorProblem(
                harm, BigFloat[1, 0], (BigFloat(0), BigFloat(1)); order)
            solb  = vector_solve_pade(probb; h = BigFloat("0.25"))
            zb = BigFloat("0.6")
            yb = solb(zb)
            @test abs(yb[1] - cos(zb)) < BigFloat(10)^(-16)
            @test abs(yb[2] + sin(zb)) < BigFloat(10)^(-16)
        end
    end

    @testset "VP.1.3 meromorphic vector system, smooth stretch" begin
        # y₁ = 1/(1−z), y₂ = 2/(1−z): both have a genuine SHARED pole at
        # z = 1.  Then y₁' = y₁², y₂' = y₁·y₂.  Integrate [0, 0.8] —
        # short of the pole — and assert accuracy against the closed
        # form.  (Pole-crossing path-network behaviour is later beads.)
        f = (z, y) -> [y[1]^2, y[1] * y[2]]
        order = 30

        prob = VectorPadeTaylorProblem(f, [1.0, 2.0], (0.0, 0.8); order)
        sol  = vector_solve_pade(prob; h = 0.2)

        for k in 1:length(sol.z)
            z = sol.z[k]
            @test sol.y[k][1] ≈ 1 / (1 - z) atol = 1e-9
            @test sol.y[k][2] ≈ 2 / (1 - z) atol = 1e-9
        end
        # Dense callable inside the meromorphic stretch.
        for z in (0.13, 0.55, 0.77)
            yz = sol(z)
            @test yz[1] ≈ 1 / (1 - z) atol = 1e-9
            @test yz[2] ≈ 2 / (1 - z) atol = 1e-9
        end
    end

    @testset "VP.1.4 dense callable correctness" begin
        # Harmonic system again; check sol(z) between breakpoints, at
        # the endpoints, and out-of-range throwing.
        harm = (z, y) -> [y[2], -y[1]]
        prob = VectorPadeTaylorProblem(harm, [1.0, 0.0], (0.0, 1.0);
                                       order = 24)
        sol  = vector_solve_pade(prob; h = 0.25)

        # Interior point not coincident with any breakpoint.
        yz = sol(0.42)
        @test yz[1] ≈ cos(0.42)  atol = 1e-8
        @test yz[2] ≈ -sin(0.42) atol = 1e-8

        # Endpoints return the IC / final state exactly.
        y_start = sol(0.0)
        @test y_start[1] ≈ 1.0 atol = 1e-14
        @test y_start[2] ≈ 0.0 atol = 1e-14
        y_end = sol(1.0)
        @test y_end[1] ≈ sol.y[end][1] atol = 1e-14
        @test y_end[2] ≈ sol.y[end][2] atol = 1e-14

        # Out-of-range z throws.
        @test_throws DomainError sol(-0.1)
        @test_throws DomainError sol(1.5)

        # Non-commensurate window — z_end is NOT a multiple of h, so the
        # final partial step MUST be clamped to land exactly on z_end.
        # (This exercises the `min(h, z_end − z)` clamp; without it the
        # trajectory overshoots and sol(z_end) is out of range.)
        prob2 = VectorPadeTaylorProblem(harm, [1.0, 0.0], (0.0, 0.9);
                                        order = 24)
        sol2  = vector_solve_pade(prob2; h = 0.25)   # 0.9 / 0.25 = 3.6
        @test sol2.z[end] ≈ 0.9 atol = 1e-14
        @test sol2.h[end] ≈ 0.9 - 3 * 0.25 atol = 1e-14   # clamped tail
        ye2 = sol2(0.9)
        @test ye2[1] ≈ cos(0.9)  atol = 1e-8
        @test ye2[2] ≈ -sin(0.9) atol = 1e-8
    end

    @testset "VP.1.5 failure modes throw informatively" begin
        harm = (z, y) -> [y[2], -y[1]]

        # h ≤ 0 — the driver rejects a non-positive step.
        prob = VectorPadeTaylorProblem(harm, [1.0, 0.0], (0.0, 1.0))
        @test_throws ArgumentError vector_solve_pade(prob; h = 0.0)
        @test_throws ArgumentError vector_solve_pade(prob; h = -0.2)

        # order < 2 — the problem constructor rejects it.
        @test_throws ArgumentError VectorPadeTaylorProblem(
            harm, [1.0, 0.0], (0.0, 1.0); order = 1)

        # Coincident zspan endpoints.
        @test_throws ArgumentError VectorPadeTaylorProblem(
            harm, [1.0, 0.0], (0.5, 0.5))

        # Empty y0.
        @test_throws ArgumentError VectorPadeTaylorProblem(
            harm, Float64[], (0.0, 1.0))

        # max_steps exceeded — a tiny step over a long window.
        prob2 = VectorPadeTaylorProblem(harm, [1.0, 0.0], (0.0, 1.0))
        @test_throws ErrorException vector_solve_pade(
            prob2; h = 0.01, max_steps = 5)
    end

    # -------------------------------------------------------------------------
    # VP.1.6 — mutation-proof.  Three mutations were applied to
    # `src/VectorProblems.jl`, the suite re-run, the RED count recorded,
    # and the mutation reverted.  This block documents the result; it is
    # not executable (the impl is in its restored, correct state).
    #
    #   M1 — wrong segment lookup in the dense callable: index the
    #        stored shared-Q at `min(k + 1, end)` instead of `k`.  sol(z)
    #        then evaluates the next segment's approximant at an
    #        out-of-range rescaled t.  Observed: 18 failures of 72 run
    #        (every VP.1.2/1.3/1.4 dense-callable assertion) — BIT.
    #   M2 — drop the final-step clamp `min(h, z_end − z)`, always step
    #        the full `h`, so a non-commensurate window overshoots
    #        z_end.  Observed: 2 failures of 72 run (VP.1.4's
    #        sol2.z[end] / sol2.h[end] clamp assertions) — BIT.  (This
    #        is why VP.1.4 carries the 0.9 / 0.25 = 3.6 non-commensurate
    #        window: commensurate windows never engage the clamp.)
    #   M3 — store the WRONG segment's Q: push the *previous* segment's
    #        denominator instead of the one just produced.  Observed:
    #        20 failures of 72 run (VP.1.2/1.3/1.4 dense callable all go
    #        off — the stored Q no longer matches its numerators) — BIT.
    #
    # All three mutations bit (counts above); impl restored to the
    # correct state and the suite re-confirmed GREEN at 72/72.
    # -------------------------------------------------------------------------

end

# =============================================================================
# VP.2.* — DESCENDING-span support (bug padetaylor-x0p0, the vector twin of the
# scalar padetaylor-xhjw fix in `Problems.solve_pade`).
#
# THE BUG.  Before the fix `vector_solve_pade` had a forward-only driver loop
# (`while real(state.z) < real(z_end)`) that NEVER ENTERED for a descending span
# (`z_end < z_start`): it silently returned a degenerate one-node trajectory
# with no error (Rule 1 violation).  The clamp (`min(h_T, z_end − state.z)`),
# the callable window guard, and the callable segment scan all baked in an
# ascending span too.
#
# THE FIX (src/VectorProblems.jl).  dir = sign(real(z_end − z_start)); the loop
# is `dir*(real(z_end) − real(state.z)) > 0`; the clamp is over POSITIVE
# magnitudes (`min(h_jz, h_T, |gap|)`) signed once by `dir`; the callable guards
# the `[lo,hi] = minmax(real(z[1]),real(z[end]))` envelope and scans with
# `cdir = sign(real(z[end]) − real(z[1]))`.  For an ascending span dir=cdir=+1
# and every expression reduces exactly to the former code (VP.1.* unchanged).
#
# ORACLE (Rule 5 — known-correct values, NOT "didn't throw").  The FW
# equianharmonic Weierstrass companion u''=6u² ⇒ y'=(y₂, 6y₁²),
# u(z)=℘(z-1; g₂=0, g₃=2), integrated LEFTWARD from z=0.  Closed-form ℘ values
# regenerated INDEPENDENTLY with wolframscript (no code shared with src/):
#   wolframscript -code 'N[WeierstrassP[#-1,{0,2}],17]&      /@ {0,-0.3,-0.5,-0.7,-1.8}'
#   wolframscript -code 'N[WeierstrassPPrime[#-1,{0,2}],17]& /@ {0,-0.3,-0.5,-1.8}'
# Nearest negative-side pole at z = 1−2ω = −1.7260681808557806
# (ω = Γ(1/3)³/(2^{13/6}π) = 1.3630340904278903): [−0.95,0] is pole-free and
# [0,−1.95] brackets exactly one pole (the descending pole-bridge edge).
# Convention check: ℘(−1;{0,2}) == u_0_FW = 1.071822516416917 (matches the FW IC).
# =============================================================================
@testset "VectorProblems descending span (x0p0)" begin
    fW_vec = (z, y) -> [y[2], 6 * y[1]^2]
    u0_FW,  up0_FW = 1.071822516416917, 1.710337353176786
    u_m0_3, up_m0_3 =  0.8012333414468928,  0.23976378516204436
    u_m0_5, up_m0_5 =  0.829689911360649,  -0.5334655722915863
    u_m0_7          =  1.0295165057128675
    u_m1_8, up_m1_8 =  182.95202472221675,  4949.209161573186

    @testset "VP.2.1 single descending segment is NOT degenerate" begin
        prob = VectorPadeTaylorProblem(fW_vec, [u0_FW, up0_FW], (0.0, -0.5);
                                       order = 30)
        sol  = vector_solve_pade(prob; h = 0.5)
        # The forward-only bug returned ZERO steps (degenerate); the fix takes one.
        @test length(sol.h) == 1
        @test sol.z[1] == 0.0
        @test sol.z[end] == -0.5            # lands on z_end exactly (clamp)
        @test sol.h[1] ≈ -0.5               # SIGNED step (dir = −1)
        # Dense callable vs closed-form ℘ (y₁ = u).
        @test isapprox(sol(0.0)[1],  u0_FW;  atol = 1e-13)
        @test isapprox(sol(-0.3)[1], u_m0_3; rtol = 1e-10)
        @test isapprox(sol(-0.5)[1], u_m0_5; rtol = 1e-10)
        # y₂ = u': a single steep segment places the derivative looser than the
        # value (cf. VPO.1 envelope); assert with a generous rtol, not ε.
        @test isapprox(sol(-0.5)[2], up_m0_5; rtol = 1e-2)
    end

    @testset "VP.2.2 multi-segment descending + pole-bridge scan" begin
        # Five clamped |h|=0.2 steps over [0,−1.0]: proves the signed clamp does
        # NOT take one giant step (the MT-2 locus).
        prob = VectorPadeTaylorProblem(fW_vec, [u0_FW, up0_FW], (0.0, -1.0);
                                       order = 30)
        sol  = vector_solve_pade(prob; h = 0.2)
        @test length(sol.h) == 5
        @test issorted(sol.z; rev = true)            # decreasing breakpoints
        @test sol.z[end] ≈ -1.0 atol = 1e-13
        @test all(h -> isapprox(h, -0.2; atol = 1e-13), sol.h)
        @test isapprox(sol(-0.7)[1], u_m0_7; rtol = 1e-10)
        @test isapprox(sol(-0.3)[1], u_m0_3; rtol = 1e-10)

        # POLE-BRIDGE: span [0,−1.95] with h=1.0 → 2 segments; segment 2
        # ([−1.0,−1.95]) brackets the pole at −1.726.  Evaluating PAST the pole
        # is load-bearing on the dir-aware SCAN (the MT-3 locus) — only segment 2
        # values z=−1.8 correctly.
        probB = VectorPadeTaylorProblem(fW_vec, [u0_FW, up0_FW], (0.0, -1.95);
                                        order = 30)
        solB  = vector_solve_pade(probB; h = 1.0)
        @test length(solB.h) == 2
        @test solB.z[end] ≈ -1.95 atol = 1e-13
        # Past the −1.726 pole.  The vector shared-Q (m,m) first-order approximant
        # bridges to ~5e-5 here — looser than the scalar 2nd-order solve_pade's
        # 1e-7 at the SAME config (problems_test.jl:164), because one shared
        # denominator must fit BOTH companion components (the documented shared-Q
        # fit characteristic, cf. the VP.1.2 entire-system residual note); NOT a
        # bug.  The assertion is SCAN-sensitive: the dir-unaware (MT-3) segment is
        # off by >0.06 here, far beyond this rtol, so it still pins the scan.
        @test isapprox(solB(-1.8)[1], u_m1_8; rtol = 1e-4)   # past the −1.726 pole
        @test isapprox(solB(-0.5)[1], u_m0_5; rtol = 1e-9)   # before the pole
    end

    @testset "VP.2.3 callable window guard (descending envelope)" begin
        prob = VectorPadeTaylorProblem(fW_vec, [u0_FW, up0_FW], (0.0, -0.5);
                                       order = 30)
        sol  = vector_solve_pade(prob; h = 0.5)
        # [lo,hi] = (−0.5, 0.0): just outside either end throws (the MT-4 locus).
        @test_throws DomainError sol(1e-9)            # above hi = 0.0
        @test_throws DomainError sol(-0.5 - 1e-9)     # below lo = −0.5
        # Just inside either end is in-window (no throw).
        @test sol(-1e-12) isa Vector
        @test sol(-0.5 + 1e-12) isa Vector
    end

    @testset "VP.2.4 forward∘reverse round-trip recovers the IC" begin
        # Self-checking reversibility invariant (no external oracle): integrate
        # ascending [0,0.5], then descending [0.5,0] from the forward endpoint;
        # the reverse solve must return to the IC.
        probF = VectorPadeTaylorProblem(fW_vec, [u0_FW, up0_FW], (0.0, 0.5);
                                        order = 30)
        solF  = vector_solve_pade(probF; h = 0.1)
        yend  = solF(0.5)
        probR = VectorPadeTaylorProblem(fW_vec, yend, (0.5, 0.0); order = 30)
        solR  = vector_solve_pade(probR; h = 0.1)
        # The load-bearing invariant is IC recovery + an EXACT landing on z_end.
        # (Segment count is float-commensurability-dependent for h=0.1: the
        # descending sum 0.5−0.1−0.1−0.1−0.1 = 0.10000000000000003 makes the
        # reverse penultimate gap exceed h by 1 ulp, so the clamp adds a final
        # sub-ε landing step — benign, and identical to the scalar driver's
        # arithmetic.  VP.2.2 pins the count on a commensurate span instead.)
        @test issorted(solR.z; rev = true)
        @test solR.z[end] == 0.0                       # final clamp lands exactly
        @test isapprox(solR(0.0)[1], u0_FW;  rtol = 1e-9)
        @test isapprox(solR(0.0)[2], up0_FW; rtol = 1e-8)
    end

    @testset "VP.2.5 BigFloat-256 descending round-trip (generic-T dir/min/minmax)" begin
        # Precision-clean: sign()/dir·min/minmax are all generic over T.  Tested
        # via the BF reversibility invariant (no BF ℘ oracle needed).
        setprecision(BigFloat, 256) do
            u0  = BigFloat("1.0718225164169174")
            up0 = BigFloat("1.7103373531767862")
            probF = VectorPadeTaylorProblem(fW_vec, BigFloat[u0, up0],
                                            (BigFloat(0), BigFloat("0.4")); order = 30)
            solF  = vector_solve_pade(probF; h = BigFloat("0.2"))
            yend  = solF(BigFloat("0.4"))
            probR = VectorPadeTaylorProblem(fW_vec, yend,
                                            (BigFloat("0.4"), BigFloat(0)); order = 30)
            solR  = vector_solve_pade(probR; h = BigFloat("0.2"))
            @test length(solR.h) == 2
            @test solR.z[end] == BigFloat(0)
            @test abs(solR(BigFloat(0))[1] - u0)  < BigFloat(10)^(-20)
            @test abs(solR(BigFloat(0))[2] - up0) < BigFloat(10)^(-18)
        end
    end

    @testset "VP.2.6 :jorba_zou adaptive descending" begin
        prob = VectorPadeTaylorProblem(fW_vec, [u0_FW, up0_FW], (0.0, -0.9);
                                       order = 30)
        sol  = vector_solve_pade(prob; h = 0.5, step_policy = :jorba_zou)
        @test length(sol.h) >= 1
        @test sol.z[end] ≈ -0.9 atol = 1e-12
        @test all(h -> real(h) < 0, sol.h)             # every step leftward
        @test all(h -> abs(h) <= 0.5 + 1e-12, sol.h)   # the h ceiling is respected
        @test isapprox(sol(-0.5)[1], u_m0_5; rtol = 1e-9)
    end
end

# -----------------------------------------------------------------------------
# VP.2.7 — mutation-proof of the x0p0 dir-aware fix (Rule 4).  Four mutations
# were applied to src/VectorProblems.jl, this file re-run in the campaign env,
# the RED count recorded, and each mutation reverted byte-clean (`git diff src/`
# empty).  Mirrors the scalar xhjw battery (MT-1..MT-4).  Results recorded after
# the serial run; this block is documentation (the impl is in its correct state).
#
# Results (campaign env /tmp/ptcampaign, 2026-06-13; VP.1.* stayed GREEN 72/72
# throughout — the ascending path is untouched):
#   MT-1 LOOP   — revert the driver loop to `while real(state.z) < real(z_end)`
#                 ⇒ descending solves take ZERO steps (degenerate one-node
#                 trajectory); 27/35 descending RED (11 fail + 16 error: the
#                 lone-node sol(z) throws DomainError as lo==hi).  BIT.
#   MT-2 CLAMP  — drop `dir *` on the step (`h_step = h_mag`, unsigned) ⇒ a
#                 descending solve steps the WRONG way and never reaches z_end:
#                 all 6 descending testsets ERROR (max_steps).  BIT.
#   MT-3 SCAN   — revert the callable scan to the ascending
#                 `real(z_T) > real(sol.z[k+1])` ⇒ 7 RED: every DESCENDING
#                 MULTI-segment dense eval picks the wrong segment — VP.2.2
#                 (sol(−0.7), sol(−0.3) AND the pole-bridge sol(−1.8)), VP.2.4,
#                 VP.2.5.  (VP.2.1 single-segment + VP.2.6's landed evals are
#                 scan-insensitive.)  BIT.
#   MT-4 GUARD  — revert the window guard to the raw `lo,hi = sol.z[1],sol.z[end]`
#                 ⇒ for a descending sol.z (lo=0 > hi) every interior negative z
#                 is rejected: 15 descending evals ERROR (DomainError).  BIT.
# All four mutations reverted byte-clean (`git diff src/VectorProblems.jl` empty).
# -----------------------------------------------------------------------------
