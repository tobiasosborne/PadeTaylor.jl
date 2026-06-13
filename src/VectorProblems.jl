"""
    PadeTaylor.VectorProblems

Public driver layer of the v0.2 vector pipeline: the problem object
`VectorPadeTaylorProblem`, the integrator entry point
`vector_solve_pade`, and the dense-output object
`VectorPadeTaylorSolution`.  This module is the vector analogue of the
scalar `Problems` module — it consumes `VectorStepper`'s
`vector_pade_step_with_pade!` segment-by-segment, collecting the
per-segment shared-`Q` approximants so the returned solution can
interpolate `y(z)` at any `z` in the integration window.

## Where this fits — the four-layer vector lift

`VectorProblems` sits one layer above `VectorStepper`, exactly as the
scalar `Problems` sits above `PadeStepper` (see `src/Problems.jl`,
"Where this fits").  `VectorStepper` knows how to take *one*
Padé–Taylor step `(z, y) → (z+h, y(z+h))` of a first-order vector ODE
`y' = f(z, y)`, `y ∈ ℂ^d`, and return the local shared-`Q` approximant
(`d` numerators over one denominator) that produced it.
`VectorProblems` owns the trajectory: it iterates the stepper,
accumulates `(z[k], y[k], h[k], shared-Q[k])` per segment, and exposes a
callable solution object.  No knowledge of Taylor jets, SVDs, or the
`h^k` rescaling lives here — those are the lower layers' concern.

## Why first-order only

The v0.2 vector pipeline integrates first-order systems `y' = f(z, y)`
*by construction*.  The scalar `Problems`/`solve_pade` carries a
1st/2nd-order branch because the v0.1 Painlevé equations are naturally
second order (`u'' = f(z, u, u')`).  v0.2's vector targets — the
Noumi–Yamada `A_n⁽¹⁾` systems and the `P_I⁽²⁾` 4-vector companion form
— are *already* first-order vector systems; any higher-order scalar ODE
is reduced to a first-order vector system before it reaches this driver.
So `VectorPadeTaylorProblem` carries no order-of-ODE branch: the state
is always the value vector `y`, never a derivative (the same essential
that makes `VectorStepper`'s state minimal — see `src/VectorStepper.jl`,
"First-order, so the state is just `(z, y)`").

## The per-segment shared-`Q` store

Each integration step produces *one* shared-denominator approximant:
`d` numerator polynomials `P₁,…,P_d` over a *single* denominator `Q`,
all in the rescaled variable `t = (z − z_seg)/h_seg ∈ [0, 1]`.  We carry
this as a small immutable `SharedPadeApproximant{T}` struct holding
`(numerators::Vector{Vector{T}}, denominator::Vector{T})` — one per
segment.  Storing the *shared* `Q` (rather than `d` independent scalar
Padés) is the v0.2 keystone: the components of a vector meromorphic
solution blow up at the *same* movable poles, so a single consistent `Q`
per segment is the honest representation (ADR-0019; `src/SharedPade.jl`,
"Why one shared denominator").  A downstream vector pole-extraction bead
(V7) reads the system's poles directly from each segment's `denominator`.

## Dense interpolation: the rescaled-variable convention

Each stored segment `k` covers `[z[k], z[k+1]]` with segment length
`h[k] = z[k+1] − z[k]`, and its shared-`Q` approximant lives in the
rescaled variable `t = (z − z[k]) / h[k] ∈ [0, 1]` — the *same* `t` the
stepper Padé-converts in (the `h^k` rescaling absorbs the step length
into the variable; `src/VectorStepper.jl`, "Why we rescale by `h^k`").
To evaluate the dense solution at an arbitrary `z` in segment `k` we set
`t = (z − z[k]) / h[k]` and return the `d`-vector

    yᵢ(z) = Pᵢ(t) / Q(t)        for component i = 1,…,d,

with every component divided by the *same* `Q(t)` — the shared-`Q`
keystone applied at evaluation time.  Unlike the scalar `Problems`
callable (which also returns `u'` via a chain-rule `1/h` factor on the
Padé derivative), the first-order vector driver returns only the *value*
`y(z)`: a first-order system carries no derivative in its state, so
there is no derivative to interpolate.

## Fail-fast contract (CLAUDE.md Rule 1)

Bad inputs throw with a `Suggestion` line: `order < 2`, coincident
`zspan` endpoints, empty `y0` (all at problem construction, mirroring
`PadeTaylorProblem`); non-positive `h` and exhausting `max_steps` (in
the driver); evaluating outside the `[lo, hi]` span envelope — for
either span ordering (in the callable; bug `padetaylor-x0p0`).
Numerical breakdowns inside the stepper — a singular shared Toeplitz
block, `Q(1) ≈ 0` meaning the step landed on a pole — propagate from
`VectorStepper`/`SharedPade` unchanged.

## Additive — distinct from the v0.1 `solve_pade`

`vector_solve_pade` is a *separate* function, NOT a method added to the
v0.1 `solve_pade` generic.  The v0.2 vector lift is additive (v0.2 plan,
"Architecture decision — additive"): the v0.1 scalar `solve_pade`
dispatch is untouched, so the v0.1 test corpus stays bit-identical.

## References

  - `src/Problems.jl` — the scalar `solve_pade` driver this module
    mirrors structurally (problem object, fixed-step loop, dense
    callable); the `d = 1` correctness oracle (VP.1.1).
  - `src/VectorStepper.jl` (commit `240a1d0`) —
    `vector_pade_step_with_pade!`, the per-step primitive this driver
    loops; the rescaled-variable convention.
  - `src/SharedPade.jl`, ADR-0019 — the shared-`Q` approximant.
  - `docs/v0p2_plan.md` — §"The 10 phases", V3b row.
"""
module VectorProblems

using ..VectorStepper:     VectorPadeStepperState, vector_pade_step_with_pade!
using ..VectorCoefficients: vector_taylor_coefficients
using ..VectorStepControl:  vector_step_jorba_zou

export VectorPadeTaylorProblem, VectorPadeTaylorSolution, vector_solve_pade

# -----------------------------------------------------------------------------
# Per-segment shared-Q approximant
# -----------------------------------------------------------------------------

"""
    SharedPadeApproximant{T}(numerators, denominator)

The shared-denominator Padé approximant for one integration segment:
`d` numerator coefficient vectors `numerators[i]` over a *single*
denominator coefficient vector `denominator`, all low-to-high and in the
rescaled variable `t ∈ [0, 1]`, normalised so `Q(0) = denominator[1]`.
Produced once per step by `VectorStepper.vector_pade_step_with_pade!`;
the dense-output callable evaluates `Pᵢ(t)/Q(t)` from it.
"""
struct SharedPadeApproximant{T}
    numerators::Vector{Vector{T}}
    denominator::Vector{T}
end

# -----------------------------------------------------------------------------
# Problem
# -----------------------------------------------------------------------------

"""
    VectorPadeTaylorProblem(f, y0::AbstractVector, zspan; order = 30)

Container for a first-order vector analytic IVP `y' = f(z, y)`,
`y(z_start) = y0`, `y ∈ ℂ^d`.  `f(z, y)` takes a scalar / `Taylor1` `z`
and a `d`-vector `y` and returns the `d`-vector `ẏ`.  The element type
`T` is promoted from the `zspan` endpoints *and* the `y0` components;
`y0` is coerced into `Vector{T}` at construction.  `order` is the Taylor
truncation degree consumed by the inner stepper (FW 2011 §5.1 settles on
30 as the canonical choice).

Throws `ArgumentError` (Rule 1) for `order < 2`, coincident `zspan`
endpoints, or an empty `y0` — mirroring the scalar `PadeTaylorProblem`'s
validation.
"""
struct VectorPadeTaylorProblem{F, T}
    f::F
    y0::Vector{T}
    zspan::Tuple{T, T}
    order::Int
end

function VectorPadeTaylorProblem(f, y0::AbstractVector, zspan::Tuple;
                                 order::Integer = 30)
    order ≥ 2 || throw(ArgumentError(
        "VectorPadeTaylorProblem: order must be ≥ 2 (got $order); the " *
        "inner Padé-Taylor stepper requires at least two Taylor passes. " *
        "Suggestion: pass `order = 30` (FW 2011 §5.1 default)."))
    isempty(y0) && throw(ArgumentError(
        "VectorPadeTaylorProblem: y0 is empty — a zero-component system " *
        "has nothing to integrate. " *
        "Suggestion: provide a non-empty initial-condition vector."))
    z_start, z_end = zspan
    z_start == z_end && throw(ArgumentError(
        "VectorPadeTaylorProblem: zspan endpoints coincide " *
        "($z_start == $z_end). " *
        "Suggestion: provide a non-degenerate interval."))
    T = promote_type(typeof(z_start), typeof(z_end), eltype(y0))
    y0_T = Vector{T}(y0)
    return VectorPadeTaylorProblem{typeof(f), T}(
        f, y0_T, (T(z_start), T(z_end)), Int(order))
end

# -----------------------------------------------------------------------------
# Solution
# -----------------------------------------------------------------------------

"""
    VectorPadeTaylorSolution{T}

Trajectory + per-segment shared-`Q` store for a vector IVP solve.  For a
trajectory of `n` segments:

  - `z::Vector{T}`             — the `n+1` breakpoints (`z[1] = z_start`,
                                 `z[end] = z_end`);
  - `y::Vector{Vector{T}}`     — the `d`-vector state `y(z[k])` at each
                                 breakpoint;
  - `h::Vector{T}`             — the `n` segment lengths
                                 `h[k] = z[k+1] − z[k]`;
  - `pade::Vector{SharedPadeApproximant{T}}` — the `n` per-segment
                                 shared-`Q` approximants, segment `k`
                                 living in `t = (z − z[k])/h[k] ∈ [0,1]`.

Callable: `sol(z) -> Vector{T}` returns the dense interpolated state
`y(z)` for any `z` in `[z_start, z_end]`; see the module docstring,
"Dense interpolation".
"""
struct VectorPadeTaylorSolution{T}
    z::Vector{T}
    y::Vector{Vector{T}}
    h::Vector{T}
    pade::Vector{SharedPadeApproximant{T}}
end

# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------

"""
    vector_solve_pade(prob::VectorPadeTaylorProblem; h, max_steps = 100_000,
                      step_policy::Symbol = :fixed)
        -> VectorPadeTaylorSolution

Take Padé–Taylor steps of the first-order vector ODE until the
integration window `zspan` is exhausted, storing each segment's
shared-`Q` approximant for later dense evaluation via the callable
interface.  The final partial step is clamped so the trajectory lands
*exactly* on `z_end`.

The span may be ASCENDING (`real(z_end) > real(z_start)`) or DESCENDING
(`real(z_end) < real(z_start)`): the driver steps with
`dir = sign(real(z_end − z_start))` and the inner stepper round-trips a
signed `h`, so a leftward integration is fully supported (bug
`padetaylor-x0p0`, the vector twin of scalar `padetaylor-xhjw`).  `h` is
always a positive-magnitude ceiling; `dir` sets the travel direction.
For a descending span `sol.z` is *decreasing*; the callable handles
either ordering.

`step_policy` selects how the step length is chosen each iteration:

  - `:fixed` (default) — every step is the supplied `h` (the final
    partial step clamped to `z_end`).  Byte-identical to the V3b
    fixed-step behaviour.
  - `:jorba_zou` — at each step the `d` Taylor jets are formed at the
    current state and `VectorStepControl.vector_step_jorba_zou` chooses
    a norm-based truncation-honouring step `h_jz`.  The supplied `h`
    acts as a **ceiling**: the actual step is
    `min(h_jz, h, z_end − z)` — never larger than the user's `h`, and
    never overshooting `z_end`.  `h` is thus the largest step the
    integrator is *allowed* to take; the adaptive policy may take a
    smaller one when the local jet demands it.

Throws `ArgumentError` for non-positive `h` or an unknown `step_policy`,
and `ErrorException` if the integration does not reach `z_end` within
`max_steps` steps (Rule 1).  Numerical breakdowns inside the stepper
propagate unchanged.
"""
function vector_solve_pade(prob::VectorPadeTaylorProblem{F, T};
                           h::Number,
                           max_steps::Integer = 100_000,
                           step_policy::Symbol = :fixed) where {F, T}
    real(h) > 0 && imag(h) == 0 || throw(ArgumentError(
        "vector_solve_pade: h must be a strictly-positive real step " *
        "length (got $h). " *
        "Suggestion: pass a positive real `h`."))
    step_policy in (:fixed, :jorba_zou) || throw(ArgumentError(
        "vector_solve_pade: unknown step_policy = :$step_policy. " *
        "Suggestion: pass `:fixed` (default) or `:jorba_zou`."))

    z_start, z_end = prob.zspan
    h_T   = T(h)
    state = VectorPadeStepperState{T}(z_start, prob.y0)
    # Integration direction (bug padetaylor-x0p0, mirroring the scalar
    # padetaylor-xhjw fix in `Problems.solve_pade`).  The inner stepper
    # round-trips a SIGNED step: `VectorStepper._rescale_by_powers` is a pure
    # `h^k` power (`src/VectorStepper.jl:291-299`) and the pole guard tests
    # `abs(Q(1))` (`:263`), both sign-agnostic, so a DESCENDING span
    # (`real(z_end) < real(z_start)`) is integrated leftward rather than
    # silently returning a degenerate one-node trajectory (Rule 1).  `dir` is
    # ±1 in `real(T)` (the constructor rejects `z_start == z_end`, so it is
    # never zero).  The user-supplied `h` stays a POSITIVE-magnitude ceiling;
    # `dir` is applied once to the clamped magnitude, so `|h_step| ≤ h` and
    # `≤ |gap|`.  For an ascending span `dir = 1` and every expression below
    # reduces exactly to the former forward-only code.
    dir = sign(real(z_end - z_start))

    z_vec    = T[z_start]
    y_vec    = Vector{T}[copy(prob.y0)]
    h_vec    = T[]
    pade_vec = SharedPadeApproximant{T}[]

    steps = 0
    while dir * (real(z_end) - real(state.z)) > 0
        steps += 1
        steps ≤ max_steps || error(
            "vector_solve_pade: did not reach z_end after " *
            "max_steps=$max_steps steps (current z=$(state.z), target " *
            "z_end=$z_end). " *
            "Suggestion: increase max_steps, or shorten the integration " *
            "window.")
        # Choose the step length.  The clamp is over POSITIVE magnitudes
        # (`h_jz` from the Jorba–Zou selector and the user's `h` are positive;
        # `gap_mag = |z_end − z|`), then signed once by `dir`.  Under :fixed it
        # is the supplied h, clamped to the gap on the final partial step so
        # the trajectory lands on `z_end` exactly.  Under :jorba_zou the
        # norm-based selector proposes `h_jz` and the supplied `h` is a
        # ceiling: `h_step = dir·min(h_jz, h, |z_end − z|)`.
        gap_mag = abs(z_end - state.z)
        if step_policy === :jorba_zou
            jets = vector_taylor_coefficients(prob.f, state.z, state.y,
                                              prob.order)
            h_jz = T(vector_step_jorba_zou(jets, eps(real(T))))
            h_mag = min(h_jz, h_T, gap_mag)
        else
            h_mag = min(h_T, gap_mag)
        end
        h_step = dir * h_mag
        _, numerators, denominator =
            vector_pade_step_with_pade!(state, prob.f, prob.order, h_step)
        push!(z_vec, state.z)
        push!(y_vec, copy(state.y))
        push!(h_vec, h_step)
        push!(pade_vec, SharedPadeApproximant{T}(numerators, denominator))
    end

    return VectorPadeTaylorSolution{T}(z_vec, y_vec, h_vec, pade_vec)
end

# -----------------------------------------------------------------------------
# Dense evaluation
# -----------------------------------------------------------------------------

"""
    (sol::VectorPadeTaylorSolution{T})(z) -> Vector{T}

Dense-output callable: return the interpolated `d`-vector state `y(z)`
for any `z` in the integration window `[z_start, z_end]`.  Locates the
segment containing `z`, rescales to `t = (z − z_seg)/h_seg`, and
evaluates each component `Pᵢ(t)/Q(t)` — every component over the *same*
`Q(t)` (the shared-`Q` keystone).  Throws `DomainError` if `z` lies
outside `[z_start, z_end]`.
"""
function (sol::VectorPadeTaylorSolution{T})(z) where {T}
    z_T = T(z)
    # Direction-agnostic window guard (bug padetaylor-x0p0): for a DESCENDING
    # span `sol.z` is decreasing, so guard against the `[lo, hi]` envelope
    # rather than the raw first/last breakpoints (which assumed ascending).
    lo, hi = minmax(real(sol.z[1]), real(sol.z[end]))
    real(z_T) < lo && throw(DomainError(z,
        "VectorPadeTaylorSolution: z=$z is outside the integration window " *
        "[$lo, $hi]."))
    real(z_T) > hi && throw(DomainError(z,
        "VectorPadeTaylorSolution: z=$z is outside the integration window " *
        "[$lo, $hi]."))

    # Linear scan for the segment containing z — v0.2 vector trajectories
    # have a modest segment count; switch to bisection if a vector
    # path-network with many segments ships later.  `cdir` makes the scan
    # advance in the direction of integration (`sol.z` may be decreasing for a
    # descending span; for an ascending span `cdir = 1` and the test reduces to
    # the former `real(z_T) > real(sol.z[k+1])`).
    cdir = sign(real(sol.z[end]) - real(sol.z[1]))
    k = 1
    @inbounds while k < length(sol.h) && cdir * (real(z_T) - real(sol.z[k + 1])) > 0
        k += 1
    end

    h_k = sol.h[k]
    t   = (z_T - sol.z[k]) / h_k
    P   = sol.pade[k]
    q_t = _eval_poly(P.denominator, t)
    return T[_eval_poly(num, t) / q_t for num in P.numerators]
end

"""
    _eval_poly(c::AbstractVector{T}, t) -> T

Horner evaluation of the polynomial `Σ c[k+1] tᵏ` (coefficients
low-to-high).  Used to evaluate each shared-`Q` numerator and the shared
denominator at the rescaled point `t`.
"""
function _eval_poly(c::AbstractVector{T}, t) where {T}
    s = zero(T)
    @inbounds for k in length(c):-1:1
        s = s * t + c[k]
    end
    return s
end

end # module VectorProblems
