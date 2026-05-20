"""
V4 end-to-end smoke test for the v0.2 vector stack (bead
`padetaylor-0ln.11`; v0.2 plan rows V2→V3a→V3b→V3c).

This is the first end-to-end validation of the whole v0.2 vector
pipeline — vector Taylor jets (`VectorCoefficients`), the shared-`Q`
stepper (`VectorStepper`), the trajectory driver (`VectorProblems`),
and the Jorba–Zou step controller (`VectorStepControl`) — against a
*genuine exact oracle*: the N=2 rational-KdV Calogero–Moser pole
system, whose pole trajectories have a closed form.

# The Calogero–Moser N=2 pole system

Krichever 1980 (`references/calogero_moser/
Krichever1980_elliptic_KP_calogero_moser_FAA14.pdf`, p. 284 eq. (12))
states the equations of motion of the integrable particle system whose
poles carry the KP/KdV elliptic solutions:

    ẍᵢ = 4 Σₖ≠ᵢ ℘'(xᵢ − xₖ),   i = 1, …, N.

The rational degeneration (Krichever p. 283 eq. (8): `℘(α) = α⁻² +
O(α²)`, hence the *repulsive* rational pole system of Airault–McKean–
Moser 1977 — poles of a real rational KdV solution repel) gives, for
N=2, the closed system

    ẍ₁ = 4/(x₁ − x₂)³,   ẍ₂ = 4/(x₂ − x₁)³.

This is the equation stated in `docs/v0p2_pillarEF_methodology_calogero
_findings.md` §Q6(b) lines 419–420.

NOTE on the sign convention.  Krichever's Hamiltonian (eq. 1, p. 282)
is *attractive* (`−2Σ℘`), and the literal `℘' = −2α⁻³` degeneration of
eq. (12) gives `ẍᵢ = −8 Σ(xᵢ−xₖ)⁻³` — an attractive system whose
symmetric IC collapses to a particle collision, with no smooth real
interval to integrate over.  The pole system relevant to *real
rational KdV solutions* (Airault–McKean–Moser 1977) is the repulsive
one, `ẍᵢ = +4/(xᵢ−xₖ)³`, stated in the pillar doc and used here: with
the symmetric IC the two poles repel along the real axis, tracing a
smooth hyperbola — exactly the smooth oracle a smoke test needs.  The
EOM constant `4` is confirmed below by independent RK4 integration.

# The verified exact oracle

Symmetric initial condition (pillar doc §Q6(b) line 532):

    x₁(0) = +r₀,  x₂(0) = −r₀,  ẋ₁(0) = ẋ₂(0) = 0.

By symmetry x₂(t) = −x₁(t) for all t, so the centre of mass stays at
0.  The relative coordinate r(t) = x₁ − x₂ = 2x₁ satisfies the single
second-order ODE

    r̈ = ẍ₁ − ẍ₂ = 4/r³ − 4/(−r)³ = 8/r³,   r(0) = 2r₀,  ṙ(0) = 0.

First integral.  `d/dt(½ṙ²) = r̈ṙ = 8ṙ/r³` and `d/dt(−4/r²) = 8ṙ/r³`,
so `½ṙ² = −4/r² + C'`, i.e. `ṙ² = C − 8/r²`.  The IC `ṙ(0)=0,
r(0)=2r₀` gives `C = 8/(2r₀)² = 2/r₀²`.  (The pillar doc §Q6(b) line
544 writes `C = 4/r₀²` and `ṙ² = C − 4/r²` — BOTH are wrong: the
coefficient is 8, not 4, and the constant follows as `2/r₀²`.  Its
final formula `x(t) = ±½√(1+4t²)` is therefore also wrong; the correct
form is derived below and is registered for a pillar-doc fix.)

Then `ṙ² = (2/r₀²)(r² − 4r₀²)/r²`, and separating variables

    r dr / √(r² − 4r₀²) = (√2/r₀) dt   ⇒   √(r² − 4r₀²) = (√2/r₀) t,

so r(t)² = 4r₀² + 2t²/r₀² and, with x₁ = r/2,

    x₁(t) = √( r₀² + t²/(2r₀²) ),   x₂(t) = −x₁(t).

The velocities follow by differentiation:

    v₁(t) = ẋ₁(t) = t / (2 r₀² x₁(t)),   v₂(t) = −v₁(t).

For r₀ = 1 this is x₁(t) = √(1 + t²/2).  This closed form was
**cross-checked independently** against a from-scratch 256-bit BigFloat
RK4 integration of the first-order CM system (200 000 steps over each
of t ∈ {0.1, 0.25, 0.4}): the RK4 trajectory agreed with √(1 + t²/2)
to 26–30 significant digits, while it *disagreed* with the pillar
doc's ±½√(1+4t²) (which gives x₁(0)=±½, contradicting the stated
IC x₁(0)=+1).  Only after that confirmation is the closed form used
below as the exact oracle.  The `CM.0` testset re-runs an abbreviated
version of that cross-check so the oracle is verified in-suite, not
trusted from a comment.

# Achievable accuracy — the honest Float64 tolerance

The v0.2 shared-`Q` stepper carries a small fixed fit residual on
smooth (entire / locally-meromorphic) systems that does not shrink
with order at Float64 precision — VP.1.2 records ~4e-12 for the
harmonic system; the CM system's rational RHS pushes it higher.
Measured worst-case Float64 errors over [0, 0.4] with order 30,
h = 0.1 are ~1.3e-10 (positions), ~6.7e-10 (velocities), ~2.7e-10
(first-integral drift).  The assertions below use `1e-8` headroom on
positions/velocities/dense output and the first integral — this is
the honest achievable accuracy of the Float64 stack, not a relaxed
tolerance to mask a bug (CLAUDE.md Rule 5 / "don't relax tolerances").
The CM.4 BigFloat run reaches ~2e-24 with the *same* code path,
confirming the residual is roundoff-limited and the stack is
genuinely precision-generic.  The centre-of-mass invariant x₁+x₂ is
*exactly* 0 — the RHS returns `[v₁, v₂, a, −a]`, so the symmetry is
structural, not numerical, and is asserted at `1e-12`.

# The conserved quantities

Two invariants are tracked along the numerical trajectory:

  - centre of mass  X(t) = x₁ + x₂ ≡ 0  (symmetric IC);
  - the relative-motion first integral  E = ½ṙ² + 4/r²  (`r̈ = 8/r³`
    is `r̈ = −∂V/∂r` with `V(r) = 4/r²`), which equals C/2 = 1/r₀² and
    is constant for the exact flow.  Equivalently in particle
    coordinates E = ½(v₁−v₂)² + 4/(x₁−x₂)².

# Mutation-proof (CLAUDE.md Rule 4; pillar doc §Q6(b) "Mutation-proof")

Two meaningful mutations were applied and the suite re-run; both bit.
Counts and the restored-GREEN confirmation are recorded in the
`CM.MUT` comment block at the end of this file.

# References

  - `references/calogero_moser/Krichever1980_elliptic_KP_calogero_
    moser_FAA14.pdf` — p. 282 eq. (1) (Hamiltonian), p. 283 eq. (8)
    (`℘ = α⁻² + O(α²)`), p. 284 eq. (12) (equations of motion).
  - `docs/v0p2_pillarEF_methodology_calogero_findings.md` §Q5, §Q6(b)
    — the CM↔KdV correspondence and the smoke-test recipe (whose
    closed-form arithmetic this file corrects).
  - `src/VectorProblems.jl`, `src/VectorStepper.jl`,
    `src/VectorStepControl.jl`, `src/VectorCoefficients.jl` — the v0.2
    vector stack under test.

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/calogero_moser_test.jl`) and under
`runtests.jl`.
"""

using Test
using PadeTaylor
using PadeTaylor.VectorProblems: VectorPadeTaylorProblem, vector_solve_pade

# -----------------------------------------------------------------------------
# The CM N=2 system as a first-order vector ODE, and the exact oracle.
# Kept test-local and independent of the impl (inline RHS — no src/ module).
# -----------------------------------------------------------------------------

# First-order state y = (x₁, x₂, v₁, v₂); ẏ = (v₁, v₂, ẍ₁, ẍ₂) with
# ẍ₁ = 4/(x₁−x₂)³, ẍ₂ = 4/(x₂−x₁)³ (repulsive rational CM, Krichever
# 1980 eq. 12 + the AMM rational-pole sign convention).  `eom_const`
# is the equation-of-motion constant (4 for the true system) — exposed
# as a parameter so the CM.MUT mutation block can perturb it.
function cm_rhs(eom_const)
    return (z, y) -> begin
        x1, x2, v1, v2 = y[1], y[2], y[3], y[4]
        d = x1 - x2
        a1 = eom_const / d^3
        return [v1, v2, a1, -a1]
    end
end

# The verified exact oracle.  r₀ generic; element type follows `t`.
cm_x1(t, r0) = sqrt(r0^2 + t^2 / (2 * r0^2))
cm_v1(t, r0) = t / (2 * r0^2 * cm_x1(t, r0))
# Full state oracle y(t) = (x₁, x₂, v₁, v₂) = (x₁, −x₁, v₁, −v₁).
cm_state(t, r0) = (cm_x1(t, r0), -cm_x1(t, r0), cm_v1(t, r0), -cm_v1(t, r0))

# Relative-motion first integral E = ½(v₁−v₂)² + 4/(x₁−x₂)².  For the
# true flow E ≡ C/2 = 1/r₀² (here r₀ = 1, so E ≡ 1).
cm_energy(y) = (y[3] - y[4])^2 / 2 + 4 / (y[1] - y[2])^2

# Independent oracle cross-check: a from-scratch RK4 integration of the
# CM first-order system.  Used ONLY in CM.0 to confirm the closed form
# before it is trusted elsewhere — never as the production oracle.
function _rk4_cm(y0::Vector{T}, t1::T, n::Int) where {T}
    rhs = cm_rhs(T(4))           # true EOM constant
    h = t1 / n
    y = copy(y0)
    for _ in 1:n
        k1 = rhs(zero(T), y)
        k2 = rhs(zero(T), y .+ (h / 2) .* k1)
        k3 = rhs(zero(T), y .+ (h / 2) .* k2)
        k4 = rhs(zero(T), y .+ h .* k3)
        y = y .+ (h / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    end
    return y
end

@testset "Calogero–Moser N=2 exact-oracle smoke test (V4)" begin

    r0   = 1.0
    f4   = cm_rhs(4.0)                 # the true CM RHS
    y0   = [r0, -r0, 0.0, 0.0]         # symmetric IC
    span = (0.0, 0.4)                  # t* (collision-free here) = ∞ for
                                       # the repulsive system; 0.4 is a
                                       # smooth interior interval.

    @testset "CM.0 closed form agrees with from-scratch RK4 (oracle check)" begin
        # Confirm the verified closed form against an INDEPENDENT
        # high-precision RK4 integration before trusting it as the
        # oracle (pillar doc §Q6(b) cross-check requirement).
        setprecision(BigFloat, 256) do
            yb0 = BigFloat[1, -1, 0, 0]
            for t in (BigFloat("0.1"), BigFloat("0.25"), BigFloat("0.4"))
                yrk = _rk4_cm(yb0, t, 120_000)
                x1e, x2e, v1e, v2e = cm_state(t, BigFloat(1))
                # RK4 with 120k steps ⇒ global error ~1e-20 or better.
                @test abs(yrk[1] - x1e) < BigFloat(10)^(-18)
                @test abs(yrk[2] - x2e) < BigFloat(10)^(-18)
                @test abs(yrk[3] - v1e) < BigFloat(10)^(-18)
                @test abs(yrk[4] - v2e) < BigFloat(10)^(-18)
            end
            # The pillar doc's formula ±½√(1+4t²) must DISAGREE — it
            # gives x₁(0) = ½, contradicting the IC x₁(0) = 1.
            bad(t) = sqrt(1 + 4 * t^2) / 2
            yrk = _rk4_cm(BigFloat[1, -1, 0, 0], BigFloat("0.25"), 120_000)
            @test abs(yrk[1] - bad(BigFloat("0.25"))) > BigFloat("0.001")
        end
    end

    @testset "CM.1 :fixed step policy matches the exact oracle" begin
        prob = VectorPadeTaylorProblem(f4, y0, span; order = 30)
        sol  = vector_solve_pade(prob; h = 0.1, step_policy = :fixed)

        # Every breakpoint: positions AND velocities match the oracle.
        for k in 1:length(sol.z)
            t = sol.z[k]
            x1e, x2e, v1e, v2e = cm_state(t, r0)
            @test sol.y[k][1] ≈ x1e atol = 1e-8
            @test sol.y[k][2] ≈ x2e atol = 1e-8
            @test sol.y[k][3] ≈ v1e atol = 1e-8
            @test sol.y[k][4] ≈ v2e atol = 1e-8
        end
        @test sol.z[end] ≈ 0.4

        # Dense callable at interior points (not on breakpoints).
        for t in (0.07, 0.13, 0.21, 0.34, 0.385)
            yt = sol(t)
            x1e, x2e, v1e, v2e = cm_state(t, r0)
            @test yt[1] ≈ x1e atol = 1e-8
            @test yt[2] ≈ x2e atol = 1e-8
            @test yt[3] ≈ v1e atol = 1e-8
            @test yt[4] ≈ v2e atol = 1e-8
        end

        # Conserved quantities along the numerical trajectory.
        for k in 1:length(sol.z)
            # Centre of mass x₁ + x₂ ≡ 0 (symmetric IC).
            @test abs(sol.y[k][1] + sol.y[k][2]) < 1e-12
            # First integral E ≡ 1/r₀² = 1.
            @test abs(cm_energy(sol.y[k]) - 1 / r0^2) < 1e-8
        end
    end

    @testset "CM.2 :jorba_zou step policy matches the exact oracle" begin
        prob = VectorPadeTaylorProblem(f4, y0, span; order = 30)
        sol  = vector_solve_pade(prob; h = 0.1, step_policy = :jorba_zou)

        for k in 1:length(sol.z)
            t = sol.z[k]
            x1e, x2e, v1e, v2e = cm_state(t, r0)
            @test sol.y[k][1] ≈ x1e atol = 1e-8
            @test sol.y[k][2] ≈ x2e atol = 1e-8
            @test sol.y[k][3] ≈ v1e atol = 1e-8
            @test sol.y[k][4] ≈ v2e atol = 1e-8
        end
        @test sol.z[end] ≈ 0.4

        for t in (0.06, 0.19, 0.28, 0.37)
            yt = sol(t)
            x1e, _, v1e, _ = cm_state(t, r0)
            @test yt[1] ≈ x1e atol = 1e-8
            @test yt[3] ≈ v1e atol = 1e-8
        end

        for k in 1:length(sol.z)
            @test abs(sol.y[k][1] + sol.y[k][2]) < 1e-12
            @test abs(cm_energy(sol.y[k]) - 1 / r0^2) < 1e-8
        end
    end

    @testset "CM.3 :fixed and :jorba_zou both reach z_end; step counts" begin
        prob   = VectorPadeTaylorProblem(f4, y0, span; order = 30)
        solf   = vector_solve_pade(prob; h = 0.1, step_policy = :fixed)
        solj   = vector_solve_pade(prob; h = 0.1, step_policy = :jorba_zou)
        nfixed = length(solf.h)
        njz    = length(solj.h)
        # :fixed with h = 0.1 over [0, 0.4] takes exactly 4 steps.
        @test nfixed == 4
        # :jorba_zou treats h = 0.1 as a ceiling; it takes ≥ 4 steps
        # (never fewer — the ceiling caps it) and reaches z_end.
        @test njz ≥ 4
        @test solf.z[end] ≈ 0.4
        @test solj.z[end] ≈ 0.4
        @info "CM.3 step counts" fixed = nfixed jorba_zou = njz
    end

    @testset "CM.4 BigFloat run reaches tighter accuracy than Float64" begin
        # The Float64 run's worst-case position error over the breakpoints.
        probf = VectorPadeTaylorProblem(f4, y0, span; order = 30)
        solf  = vector_solve_pade(probf; h = 0.1, step_policy = :fixed)
        err_f64 = maximum(
            abs(solf.y[k][1] - cm_x1(solf.z[k], r0))
            for k in 1:length(solf.z))

        # The same solve at 256-bit BigFloat precision.
        err_big = setprecision(BigFloat, 256) do
            f4b  = cm_rhs(BigFloat(4))
            y0b  = BigFloat[1, -1, 0, 0]
            spb  = (BigFloat(0), BigFloat("0.4"))
            probb = VectorPadeTaylorProblem(f4b, y0b, spb; order = 30)
            solb  = vector_solve_pade(probb; h = BigFloat("0.1"),
                                      step_policy = :fixed)
            maximum(abs(solb.y[k][1] - cm_x1(solb.z[k], BigFloat(1)))
                    for k in 1:length(solb.z))
        end

        # Float64 reaches ~1e-10; BigFloat (256-bit) reaches ~2e-24 with
        # the SAME code path — far tighter, demonstrating the
        # precision-generic stack.  The bounds carry honest headroom on
        # the measured values (CLAUDE.md Rule 5).
        @test err_f64 < 1e-9
        @test err_big < BigFloat(10)^(-22)
        @test err_big < err_f64 / 1e6
        @info "CM.4 accuracy" float64 = err_f64 bigfloat = Float64(err_big)
    end

end

# -----------------------------------------------------------------------------
# CM.MUT — mutation-proof (CLAUDE.md Rule 4; pillar doc §Q6(b)
# "Mutation-proof instruction").  Two meaningful mutations were applied,
# the suite re-run, the RED count recorded, and the mutation reverted.
# This block documents the result; it is not executable (impl + oracle
# are in their restored, correct state).
#
#   M1 — wrong EOM exponent.  In `cm_rhs`, change the acceleration
#        `eom_const / d^3` to `eom_const / d^2` (the pillar doc's
#        suggested mutation: ẍ = 4/r³ → 4/r²).  The v0.2 stepper then
#        integrates the WRONG vector field; its trajectory no longer
#        matches the √(1 + t²/2) oracle, and the first integral
#        E = ½ṙ² + 4/r² is no longer conserved.  `cm_rhs` also backs
#        the CM.0 RK4 cross-check, so that testset's RK4 trajectory is
#        likewise wrong.
#        Observed: 83 failures of 110 tests (CM.0's 12 RK4-vs-closed-
#        form comparisons, every CM.1/CM.2 position/velocity/dense/
#        energy assertion, all 3 CM.4 accuracy bounds) — BIT.  CM.3
#        (step-count only; integration still reaches z_end) stayed
#        GREEN.
#
#   M2 — wrong oracle formula.  In `cm_x1`, change
#        `sqrt(r0^2 + t^2/(2*r0^2))` to `sqrt(r0^2 + t^2/r0^2)` (drop
#        the factor ½ in the t² term — the pillar-doc-style slip).
#        The v0.2 stepper is correct but is checked against a wrong
#        oracle, so every comparison and the CM.0 RK4 cross-check fail.
#        Observed: 75 failures of 110 tests (CM.0's 12 oracle-check
#        comparisons, all CM.1/CM.2 position/velocity/dense assertions,
#        all 3 CM.4 accuracy bounds) — BIT.  The centre-of-mass and
#        first-integral assertions (which do not reference `cm_x1`) and
#        CM.3 stayed GREEN, correctly localising the mutation.
#
# Both mutations bit (counts above); `cm_rhs` and `cm_x1` were restored
# to their correct state and the suite re-confirmed GREEN at 110/110.
# -----------------------------------------------------------------------------
