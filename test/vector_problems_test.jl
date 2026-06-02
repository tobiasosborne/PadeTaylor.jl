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
