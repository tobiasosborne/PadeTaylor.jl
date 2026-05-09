"""
    PadeTaylor.Problems

Public driver layer of the four-layer architecture (ADR-0001): the
problem object `PadeTaylorProblem`, the integrator entry point
`solve_pade`, and the dense-output object `PadeTaylorSolution`.  This
module consumes `PadeStepper.pade_step_with_pade!` segment-by-segment,
collecting the per-segment Padé approximants so the returned solution
can interpolate `u` and `u'` at any `z` in the integration window.

## Where this fits — the four-layer architecture

`Problems` sits one layer above `PadeStepper`.  `PadeStepper` knows how
to take *one* Padé-Taylor step `(z, u, u') → (z+h, u(z+h), u'(z+h))`
and return the local Padé approximant `P_u` that produced it.
`Problems` owns the trajectory: it iterates the stepper, accumulates
`(z[k], y[k], h[k], pade[k])` tuples, and exposes a callable solution
object.  No knowledge of Taylor coefficients, SVDs, or rescaling lives
here — those are all the lower layers' concerns.

## The Padé-vs-Taylor pole-bridge demonstration (worklog 004)

Phase-6 v1 acceptance is a *demonstration* rather than a long-range
integration.  The canonical test problem is `u'' = 6u^2` with the
Fornberg-Weideman 2011 ICs at `z = 0`; the closed-form solution
`u(z) = ℘(z + c₁; 0, c₂)` has a lattice pole at `z = 1`.  With
`h_max = 1.5` we build a *single* segment that brackets the pole
(at rescaled `t = 1/1.5 ≈ 0.667`, strictly inside `[0, 1]`).  The same
stored Padé, evaluated at `t = z / 1.5`, returns correct values both
*before* the pole (`z = 0.5, 0.95`) and *after* (`z = 1.05, 1.4`).
A side-by-side `taylor_eval` of the same Taylor coefficients diverges
past the natural radius of convergence at `z = 1`; Padé does not.
This is the analytic-continuation advantage of the algorithm in one
test setup; tests 6.1.3 + 6.1.5 are the headline.

## v1 / v2 acceptance scope (worklog 004 + bead `padetaylor-8cr`)

**v1 (this implementation)**: pure fixed-`h_max` stepping, no vault,
no Jorba-Zou.  The architecture supports *multi*-segment trajectories
(the `solve_pade` loop iterates while `state.z < z_end`), but the v1
test corpus only exercises the single-segment `h_max = 1.5` case.

**v2 (deferred, P0 bead `padetaylor-8cr`)**: the FW 2011 §3.1 path-
network for long-range integration crossing many lattice poles, and/or
adaptive step-size selectors that compose with Padé-bridge stepping.
Two prior subagents established that fixed-`h_max` real-axis stepping
alone cannot reach the FW 2011 Table 5.1 `5e-13` budget at `z = 30`;
the path-network is the only route.  See worklog 004 for the failure
analysis.

## Dense interpolation: `t = (z - z[k]) / h[k]`,  `u' = P'(t) / h[k]`

Each stored `pade[k]` is the local Padé built in the *rescaled*
variable `t ∈ [0, 1]` covering the segment `[z[k], z[k+1]]`, with
`h[k] = z[k+1] - z[k]`.  To evaluate at an arbitrary `z` in that
segment we set `t = (z - z[k]) / h[k]`.  For `u` this is just
`_evaluate_pade(pade[k], t)`.  For `u'` the chain rule on the rescaling
contributes the `1/h[k]` factor:

    u(z[k] + h')  = P_k(h' / h[k])
    u'(z[k] + h') = P_k'(h' / h[k]) · (1 / h[k])

This is the same idiom `pade_step_with_pade!` uses to recover `u'` at
the segment endpoint (PadeStepper.jl §"Why analytic differentiation
beats re-Padé"); we generalise to arbitrary `t ∈ [0, 1]` here.

## Fail-fast contract (CLAUDE.md Rule 1)

Bad inputs throw with a `Suggestion` line: negative `order`, empty
`zspan`, non-positive `h_max`, evaluating outside `[z_start, z_end]`,
exhausting `max_steps` before reaching `z_end`.  Numerical breakdowns
inside the stepper (singular `C̃`, `Q(t) = 0` at evaluation point)
propagate from the lower layers unchanged.

## References

  - `docs/worklog/004-phase-6-pivot.md` — Phase-6 spec and oracle plan.
  - `docs/worklog/003-phase-5-padestepper.md` — `pade_step_with_pade!`
    primitive that this layer composes.
  - `src/PadeStepper.jl` — the inner-loop integration.
  - FW 2011 §3.2 line 396 — the `h^k` rescaling rationale (consumed by
    `PadeStepper`; surfaces here only via the `1/h_k` chain-rule
    factor on `u'`).
"""
module Problems

using ..RobustPade:   PadeApproximant
using ..PadeStepper:  PadeStepperState, pade_step_with_pade!,
                      _evaluate_pade, _evaluate_pade_deriv

export PadeTaylorProblem, solve_pade, PadeTaylorSolution, taylor_eval

# -----------------------------------------------------------------------------
# Problem
# -----------------------------------------------------------------------------

"""
    PadeTaylorProblem(f, y0, zspan; order = 30)

Container for an analytic IVP.  `y0 isa Tuple{T,T}` selects the
2nd-order branch (`u'' = f(z, u, u')`); a scalar `y0::T` selects the
1st-order branch (`y' = f(z, y)`).  The element type `T` is promoted
from the `zspan` endpoints; `y0` is coerced into that type at
construction.  `order` is the Taylor truncation degree consumed by the
inner stepper (FW 2011 §5.1 line 277 settles on 30 as the canonical
choice for `h ≈ 0.5` near the Painlevé-I tritronquée wall).
"""
struct PadeTaylorProblem{F, T, Y}
    f::F
    y0::Y
    zspan::Tuple{T, T}
    order::Int
end

function PadeTaylorProblem(f, y0, zspan::Tuple; order::Integer = 30)
    order ≥ 2 || throw(ArgumentError(
        "PadeTaylorProblem: order must be ≥ 2 (got $order); the inner " *
        "Padé-Taylor stepper requires at least two Taylor passes. " *
        "Suggestion: pass `order = 30` (FW 2011 §5.1 default)."))
    z_start, z_end = zspan
    T = promote_type(typeof(z_start), typeof(z_end))
    z_start == z_end && throw(ArgumentError(
        "PadeTaylorProblem: zspan endpoints coincide ($z_start == $z_end). " *
        "Suggestion: provide a non-degenerate interval."))
    if y0 isa Tuple
        y0_T = (T(y0[1]), T(y0[2]))
        return PadeTaylorProblem{typeof(f), T, typeof(y0_T)}(
            f, y0_T, (T(z_start), T(z_end)), Int(order))
    else
        y0_T = T(y0)
        return PadeTaylorProblem{typeof(f), T, typeof(y0_T)}(
            f, y0_T, (T(z_start), T(z_end)), Int(order))
    end
end

# -----------------------------------------------------------------------------
# Solution
# -----------------------------------------------------------------------------

"""
    PadeTaylorSolution{T, Y, P}

Trajectory + per-segment Padé store.  `z[k]` for `k = 1:n+1` are the
segment breakpoints (`z[1] = z_start`, `z[end] = z_end`); `y[k]` is the
state at `z[k]`; `h[k] = z[k+1] - z[k]` for `k = 1:n`; and `pade[k]`
is the local Padé approximant covering segment `k` in the rescaled
variable `t = (z - z[k]) / h[k] ∈ [0, 1]`.

Callable: `sol(z) -> (u, u')` for the 2nd-order branch, `sol(z) -> u`
for the 1st-order branch (currently unsupported in v1; see module
docstring).
"""
struct PadeTaylorSolution{T, Y, P}
    z::Vector{T}
    y::Vector{Y}
    h::Vector{T}
    pade::Vector{P}
end

# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------

"""
    solve_pade(prob::PadeTaylorProblem; h_max, max_steps = 100_000)
        -> PadeTaylorSolution

Take fixed-`h_max` Padé-Taylor steps until the integration window is
exhausted.  Each segment stores the local Padé approximant for later
dense evaluation via the callable interface.
"""
function solve_pade(prob::PadeTaylorProblem{F, T, Y};
                    h_max::Real,
                    max_steps::Integer = 100_000) where {F, T, Y}
    h_max > 0 || throw(ArgumentError(
        "solve_pade: h_max must be positive (got $h_max). " *
        "Suggestion: pass a strictly-positive step length."))
    Y <: Tuple || error(
        "solve_pade: 1st-order (scalar y0) branch is not implemented in " *
        "v1. Suggestion: rewrite the problem as a 2nd-order system and " *
        "pass `y0 = (u0, up0)`, or file a bead requesting 1st-order " *
        "support.")

    z_start, z_end = prob.zspan
    h_max_T = T(h_max)
    state   = PadeStepperState{T}(z_start, prob.y0[1], prob.y0[2])

    P_T = PadeApproximant{T}
    z_vec    = T[z_start]
    y_vec    = Y[prob.y0]
    h_vec    = T[]
    pade_vec = P_T[]

    steps = 0
    while state.z < z_end
        steps += 1
        steps ≤ max_steps || error(
            "solve_pade: did not reach z_end after max_steps=$max_steps " *
            "steps (current z=$(state.z), target z_end=$z_end). " *
            "Suggestion: increase max_steps, or shorten the integration " *
            "window.")
        h_step = min(h_max_T, z_end - state.z)
        _, P_u = pade_step_with_pade!(state, prob.f, prob.order, h_step)
        push!(z_vec, state.z)
        push!(y_vec, (state.u, state.up))
        push!(h_vec, h_step)
        push!(pade_vec, P_u)
    end

    return PadeTaylorSolution{T, Y, P_T}(z_vec, y_vec, h_vec, pade_vec)
end

# -----------------------------------------------------------------------------
# Dense evaluation
# -----------------------------------------------------------------------------

function (sol::PadeTaylorSolution{T, Y, P})(z) where {T, Y, P}
    z_T = T(z)
    z_T < sol.z[1] && throw(DomainError(z,
        "PadeTaylorSolution: z=$z is below z_start=$(sol.z[1])."))
    z_T > sol.z[end] && throw(DomainError(z,
        "PadeTaylorSolution: z=$z is above z_end=$(sol.z[end])."))

    # Why: linear scan is fine — v1 trajectories have ≤ a handful of
    # segments. Switch to bisection when the path-network ships in v2.
    k = 1
    @inbounds while k < length(sol.h) && z_T > sol.z[k + 1]
        k += 1
    end

    h_k = sol.h[k]
    t   = (z_T - sol.z[k]) / h_k
    P_u = sol.pade[k]
    u   = _evaluate_pade(P_u, t)
    up  = _evaluate_pade_deriv(P_u, t) / h_k
    return (u, up)
end

# -----------------------------------------------------------------------------
# Plain-Taylor evaluator (for the side-by-side comparison test 6.1.5)
# -----------------------------------------------------------------------------

"""
    taylor_eval(coefs::Vector{T}, h::Real) -> T

Horner-evaluate the truncated Taylor polynomial `Σ_{k=0}^{N} coefs[k+1]·h^k`
at `h`.  Used by the headline pole-bridge demo (test 6.1.5) to show
that, given the *same* Taylor coefficients, plain truncation diverges
past the natural radius of convergence at the nearest pole while the
Padé conversion does not.
"""
function taylor_eval(coefs::Vector{T}, h::Real) where {T}
    h_T = T(h)
    out = zero(T)
    @inbounds for k in length(coefs):-1:1
        out = out * h_T + coefs[k]
    end
    return out
end

end # module Problems
