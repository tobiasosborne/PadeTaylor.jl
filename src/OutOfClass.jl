"""
    PadeTaylor.OutOfClass

Enforce the *meromorphic-only* contract of `solve_pade` (CLAUDE.md Rule 1:
fail loud, never silently lie).  This module is the guard that turns the
`padetaylor-v1ub` silent-wrong-answer bug into a fail-loud
`OutOfClassError`: when `solve_pade` is driven toward a NON-pole
singularity (an essential singularity, a natural boundary — anything the
de Montessus / Nuttall–Pommerenke convergence theory does not cover), the
solver now refuses with context instead of returning a finite, plausible,
*wrong* value.

## The bug this module closes (padetaylor-v1ub)

The motivating out-of-class input is the essential-singularity recast

    u'' = u·(1 + 2z) / z⁴        exact   u(z) = e^{1/z}

whose `z = 0` is an essential singularity (the Laurent series
`1 + 1/z + 1/(2z²) + …` has infinitely many negative-power terms — NOT a
pole).  Marching `solve_pade` from `z₀ = -1` toward `0` at fixed `h = 0.1`,
the relative error against the closed-form oracle (mpmath dps=50) is

    |z|=0.5  relerr 2e-16   |z|=0.1  relerr 4e-10
    |z|=0.05 relerr 4.9e-2  |z|=0.02 relerr 1.2e17  AND WRONG SIGN

— finite throughout, no throw, no NaN.  e^{1/z} > 0 for all real z, yet
the solver returns `u ≈ -2.4e-5` near `z = 0` with no complaint.  That is
a pure Rule-1 violation in a library whose entire value proposition is
"honest results."  The empirical write-up and the auto-flip marker live in
`test/corpus_out_of_class_test.jl` (bead `padetaylor-v1ub`, ADR-0033).

## Why the discriminator is *Padé convergence*, not *radius*

The obvious-looking guard — watch the local Taylor jet's coefficient
growth, or the Jorba–Zou radius-of-convergence estimate, and trip when the
radius collapses — DOES NOT WORK, and getting this wrong would break the
package's headline feature (bridging poles).  The radius-of-convergence
equals the distance to the nearest singularity for BOTH a pole and an
essential point, so it collapses identically in both cases.  A radius
guard cannot tell "I am approaching a pole I can bridge" from "I am
approaching an essential singularity I cannot."

The correct discriminator is *Padé sequence convergence*.  GGT 2013 §8
(`references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/
GGT2013_robust_pade_via_SVD_SIREV55.md`, the "Modified Padé Approximant
with Pointwise Convergence" section, ~md:386) lays out the ground truth:

  - For a **meromorphic** function the Padé approximants converge (in
    capacity / in measure — the Nuttall–Pommerenke theorem [1,21,26]),
    so successive-order approximants AGREE.  This holds *including while
    bridging a pole*: a diagonal Padé arranges its denominator to vanish
    at the pole, and the [m,n] and [m-1,n-1] approximants of the SAME jet
    land on the same value at the step endpoint.  This convergence is
    exactly why the package can see past a pole.

  - For an **essential** singularity (GGT 2013 Fig. 9 is precisely
    `exp((z+1.5)⁻²)`, an essential point) NO finite-degree rational
    converges — the Padé machinery is not even ill-posed-but-stable, it is
    simply outside the convergence theory.  Successive-order approximants
    DISAGREE, and the disagreement grows as the jet senses the
    essential singularity more strongly.

So the per-step discriminator is the **two-order Padé convergence
defect**: from the SAME scaled Taylor jet the step already builds, form
the Padé at the working order `[m, n]` (m = n = order ÷ 2) AND at the
reduced order `[m-1, n-1]` (truncate the scaled jet by two coefficients —
NO extra ODE/jet evaluation), evaluate BOTH at the step endpoint `t = 1`,
and take the scale-relative disagreement

    δ = |P(1) − P_lo(1)| / (|P(1)| + tiny).

For an in-class (meromorphic) step δ sits at the rational-approximation
floor (machine precision while marching; ≤ ~4e-7 even for the hardest
single-step double-pole bridge).  For the e^{1/z} approach δ climbs
`…2.7e-14 → 3.2e-10 → 0.99` over the final steps.  The two populations
separate by ~6 orders of magnitude (ADR-0033 separation table).

## Why a *history gate*, not a single-step threshold

A single-step `δ > τ` test would be unsafe: the headline "bridge across
z = 0 in one big step" (CFail.1e) and the single-segment pole bridges
(CFail / CPB) produce an *isolated* δ that, while still below τ here, must
never be allowed to trip the guard on a one-off transient.  The guard
therefore fires only on **sustained, monotone growth**: δ must exceed the
threshold τ AND have grown monotonically across the last `K` step-to-step
transitions (K ≥ 2).  This is what makes the guard safe:

  - It FIRES on the e^{1/z} approach — δ grows monotonically over the last
    several steps (`2.7e-14 < 3.2e-10 < 0.99`), the signature of a jet
    losing meromorphy as it nears the essential singularity.

  - It CANNOT fire on the single-step across-0 bridge (CFail.1e): one step
    leaves no ≥2-step history, so the monotone-growth gate is never
    satisfied.

  - It IGNORES transient single-step δ spikes (an isolated pole-bridge
    step): no sustained history ⇒ no fire, and in any case those spikes
    sit far below τ.

## Calibrated thresholds (ADR-0033 — DERIVED from data, not fitted)

  - `τ = 1e-3`.  Max measured in-class δ (across every marching and
    single-step pole-bridge case in the regression corpus) is ~4e-7;
    the e^{1/z} guard must fire at δ ≈ 0.99.  τ = 1e-3 sits in that gap
    with ~3.8 orders of margin above the in-class ceiling and ~3 orders
    below the firing value.

  - `K = 2`.  The e^{1/z} approach gives ≥2 consecutive monotone-growth
    transitions ending at the firing step; the across-0 single step and
    isolated single-segment bridges give none.

## Default-on, with an escape hatch

`solve_pade(...; check_in_class = true)` enables the guard by default
(Rule 1: fail loud).  A user who is *knowingly* probing an out-of-class
input — or who has an input they are certain is meromorphic but that
trips a false positive — can pass `check_in_class = false` to restore the
legacy unguarded behaviour.  The escape hatch is documented on
`solve_pade`.

## References

  - GGT 2013 §8 (Nuttall–Pommerenke / de Montessus de Ballore; Fig. 9
    essential-singularity example) —
    `references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/
    GGT2013_robust_pade_via_SVD_SIREV55.md` (~md:386).
  - ADR-0033 — the decision, the separation table, the τ/K calibration.
  - bug `padetaylor-v1ub`; `test/corpus_out_of_class_test.jl` (corpus).
  - `src/PadeStepper.jl` — the step primitives reused here
    (`_rescale_by_powers`, `_evaluate_pade`); the unchecked
    `pade_step_with_pade!` is left byte-for-byte unchanged so the
    perf-critical `path_network_solve` hot path pays nothing for the
    guard, and `OutOfClass` provides a SEPARATE checked stepper.
  - `src/VectorWedgeStep.jl:211` — `VectorWalkError`, the typed-exception
    idiom this module's `OutOfClassError` follows.
"""
module OutOfClass

using ..Coefficients: taylor_coefficients_2nd
using ..RobustPade:   PadeApproximant, robust_pade
using ..PadeStepper:  PadeStepperState, _rescale_by_powers, _evaluate_pade

export OutOfClassError, OutOfClassChecker,
       pade_step_with_defect!, check_in_class!,
       OUT_OF_CLASS_TAU, OUT_OF_CLASS_K

# -----------------------------------------------------------------------------
# Calibrated thresholds (ADR-0033; DERIVED from the separation probe).
# -----------------------------------------------------------------------------

"""Convergence-defect threshold τ (ADR-0033): δ must exceed this to fire.
Max in-class δ ≈ 4e-7; e^{1/z} fires at δ ≈ 0.99; τ = 1e-3 sits in the gap."""
const OUT_OF_CLASS_TAU = 1.0e-3

"""Monotone-growth gate K (ADR-0033): δ must have grown across ≥ K consecutive
step-to-step transitions.  K = 2 fires on the e^{1/z} approach, never on a
single isolated bridge step."""
const OUT_OF_CLASS_K = 2

# `(|P(1)| + tiny)` floor: keeps δ finite if the endpoint value underflows to
# exactly zero.  Far below any meaningful magnitude, so it never perturbs δ.
const _DEFECT_TINY = 1.0e-300

# -----------------------------------------------------------------------------
# Typed error
# -----------------------------------------------------------------------------

"""
    OutOfClassError <: Exception

Thrown by `solve_pade` when the integrand appears to leave the documented
*meromorphic-only* class — i.e. the two-order Padé convergence defect δ has
grown monotonically and crossed the calibrated threshold τ, the signature
of a non-pole (essential) singularity for which no finite-degree rational
converges (GGT 2013 §8; see the module chapter).

Fields:

  - `z`     — the independent-variable location at which the guard tripped.
  - `δ`     — the convergence defect at that step (≫ τ).
  - `δhist` — the recent δ history (the monotone-growth run that fired).
  - `msg`   — the full human-readable message, retaining the `Suggestion:`
    remediation line (Rule 1: a fail-loud throw carries actionable advice).

`Base.showerror` prints `msg` verbatim.  The dedicated type is
self-documenting and `catch`-able (`catch e; e isa OutOfClassError`),
mirroring the `VectorWalkError` idiom (`src/VectorWedgeStep.jl:211`).
"""
struct OutOfClassError <: Exception
    z     :: Any
    δ     :: Float64
    δhist :: Vector{Float64}
    msg   :: String
end

Base.showerror(io::IO, e::OutOfClassError) = print(io, e.msg)

function _out_of_class_error(z, δ, δhist)
    trend = join((string(round(d; sigdigits = 3)) for d in δhist), " → ")
    msg = string(
        "OutOfClassError: the ODE appears to be OUT OF CLASS at z=$z.  ",
        "solve_pade is MEROMORPHIC-ONLY (it bridges poles, not essential ",
        "singularities or natural boundaries); the two-order Padé ",
        "convergence defect δ has grown monotonically and reached ",
        "δ=$(round(δ; sigdigits = 4)) (threshold τ=$(OUT_OF_CLASS_TAU)).  ",
        "For a meromorphic solution the Padé sequence converges (GGT 2013 ",
        "§8, Nuttall–Pommerenke), so δ stays near the floor even while ",
        "bridging a pole; sustained monotone growth of δ is the signature ",
        "of a non-pole singularity for which no finite-degree rational ",
        "converges.  Recent δ trend: $trend.  ",
        "Suggestion: if this input is genuinely meromorphic, pass ",
        "`check_in_class = false` to disable the guard; otherwise this ODE ",
        "is outside the documented class — recast it (e.g. an exponential ",
        "coordinate map onto a domain free of the essential singularity) ",
        "or restrict the integration window to stay clear of z=$z.")
    return OutOfClassError(z, δ, copy(δhist), msg)
end

# -----------------------------------------------------------------------------
# Two-order Padé convergence defect from a scaled jet
# -----------------------------------------------------------------------------

"""
    two_order_defect(scaled::Vector{T}, m::Int, n::Int) -> Real

The per-step in-class convergence defect δ (see the module chapter).  Given
the SCALED Taylor jet `c̃_k = h^k c_k` (the exact vector the step builds)
and the working diagonal order `m = n = order ÷ 2`, build the working-order
Padé `P = robust_pade(scaled, m, n)` AND the reduced-order
`P_lo = robust_pade(scaled[1:end-2], m-1, n-1)`, evaluate both at the step
endpoint `t = 1`, and return

    δ = |P(1) − P_lo(1)| / (|P(1)| + tiny).

No extra ODE/jet evaluation — the reduced Padé reuses the same jet,
truncated by two coefficients.  Returns a real-typed magnitude (the input
may be complex in path-network mode).  When `m ≤ 0` (degenerate tiny order)
the reduced order is clamped at `[0, 0]` so the call is always well-defined.
"""
function two_order_defect(scaled::AbstractVector{T}, m::Int, n::Int) where {T}
    one_T = one(T)
    P  = robust_pade(scaled, m, n)
    P1 = _evaluate_pade(P, one_T)
    m_lo = max(m - 1, 0)
    n_lo = max(n - 1, 0)
    # Truncate the scaled jet by two coefficients for the reduced Padé; guard
    # the tiny-order edge so the slice always has ≥ m_lo + n_lo + 1 entries.
    lo_len  = min(length(scaled) - 2, length(scaled))
    lo_len  = max(lo_len, m_lo + n_lo + 1)
    lo_len  = min(lo_len, length(scaled))
    P_lo  = robust_pade(@view(scaled[1:lo_len]), m_lo, n_lo)
    P1_lo = _evaluate_pade(P_lo, one_T)
    return abs(P1 - P1_lo) / (abs(P1) + _DEFECT_TINY)
end

# -----------------------------------------------------------------------------
# Checked stepper: one step + its convergence defect
# -----------------------------------------------------------------------------

"""
    pade_step_with_defect!(state, f, order, h) -> (state, P_u, δ::Real)

Take one Padé-Taylor step (identical numerics to
`PadeStepper.pade_step_with_pade!`) AND return the two-order convergence
defect δ for that step.  This is the SEPARATE checked entry point the
driver uses when `check_in_class = true`; the unchecked
`pade_step_with_pade!` is left byte-for-byte unchanged so the perf-critical
`path_network_solve` hot path pays nothing for the guard (and the
`@inferred`/allocation type-stability gate on the unchecked stepper is
untouched).

The step body mirrors `pade_step_with_pade!` (build the order-`order`
Taylor jet, rescale by `h^k`, diagonal `[m, n]` Padé with `m = n =
order ÷ 2`, evaluate `u` and `u'` at `t = 1`, mutate state).  The only
addition is the reduced-order `[m-1, n-1]` Padé and the δ it yields — one
extra `robust_pade` plus two Horner evals on the jet already in hand, no
extra ODE/jet evaluation.
"""
function pade_step_with_defect!(state::PadeStepperState{T}, f,
                                order::Int, h::Number) where {T}
    order ≥ 2 || throw(ArgumentError(
        "pade_step_with_defect!: order must be ≥ 2 (got $order)."))
    h_T = T(h)
    coefs_u        = taylor_coefficients_2nd(f, state.z, state.u, state.up, order)
    coefs_u_scaled = _rescale_by_powers(coefs_u, h_T)
    m = n = order ÷ 2
    P_u = robust_pade(coefs_u_scaled, m, n)
    δ   = two_order_defect(coefs_u_scaled, m, n)
    one_T  = one(T)
    new_u  = _evaluate_pade(P_u, one_T)
    new_up = _evaluate_pade_deriv_via(P_u, one_T) / h_T
    state.z  = state.z + h_T
    state.u  = new_u
    state.up = new_up
    return state, P_u, δ
end

# The checked stepper recomputes u' by the same quotient-rule sweep as
# PadeStepper._evaluate_pade_deriv.  We re-derive it locally rather than
# import the private symbol, keeping the checked path self-contained.
function _evaluate_pade_deriv_via(P::PadeApproximant{T}, z::T) where {T}
    n = length(P.a); m = length(P.b)
    N = zero(T);  @inbounds for k in n:-1:1; N = N * z + P.a[k]; end
    D = zero(T);  @inbounds for k in m:-1:1; D = D * z + P.b[k]; end
    Nt = zero(T); @inbounds for k in n:-1:2; Nt = Nt * z + (k - 1) * P.a[k]; end
    Dt = zero(T); @inbounds for k in m:-1:2; Dt = Dt * z + (k - 1) * P.b[k]; end
    iszero(D) && throw(DomainError(z,
        "pade_step_with_defect!: local Padé denominator vanishes at the step " *
        "endpoint.  Suggestion: shorten the step length."))
    return (Nt * D - N * Dt) / (D * D)
end

# -----------------------------------------------------------------------------
# History-gated checker
# -----------------------------------------------------------------------------

"""
    OutOfClassChecker(; τ = OUT_OF_CLASS_TAU, K = OUT_OF_CLASS_K)

Stateful per-solve guard.  Holds the recent δ history (a short ring buffer
of length `K + 1`, enough to test `K` consecutive monotone-growth
transitions) and the threshold/gate parameters.  Constructed once before
the `solve_pade` loop; fed each step's δ via `check_in_class!`.
"""
mutable struct OutOfClassChecker
    τ       :: Float64
    K       :: Int
    history :: Vector{Float64}   # most-recent-last, capped at length K+1
end

OutOfClassChecker(; τ::Real = OUT_OF_CLASS_TAU, K::Integer = OUT_OF_CLASS_K) =
    OutOfClassChecker(Float64(τ), Int(K), Float64[])

"""
    check_in_class!(checker, δ, z)

Record this step's defect δ and apply the throw policy.  Throws
`OutOfClassError` iff `δ > checker.τ` AND the last `K + 1` recorded δ
values are strictly increasing (i.e. δ has grown monotonically across `K`
consecutive step-to-step transitions).  Otherwise returns `nothing` and the
solve continues.

The monotone-growth gate is the safety mechanism (see the module chapter):
a single isolated δ spike — the across-0 bridge, a one-segment pole bridge
— leaves no ≥K-step history and therefore cannot trip the guard, while the
sustained climb of the e^{1/z} approach does.
"""
function check_in_class!(checker::OutOfClassChecker, δ::Real, z)
    δf = Float64(δ)
    push!(checker.history, δf)
    cap = checker.K + 1
    length(checker.history) > cap && popfirst!(checker.history)

    δf > checker.τ || return nothing
    length(checker.history) ≥ cap || return nothing   # not enough history yet
    @inbounds for i in 2:length(checker.history)
        checker.history[i] > checker.history[i - 1] || return nothing
    end
    throw(_out_of_class_error(z, δf, checker.history))
end

end # module OutOfClass
