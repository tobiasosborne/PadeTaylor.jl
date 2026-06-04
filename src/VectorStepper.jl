"""
    PadeTaylor.VectorStepper

Take one Padé–Taylor step of a first-order **vector** analytic ODE
`y' = f(z, y)`, `y ∈ ℂ^d`, advancing the state `(z, y)` to
`(z+h, y(z+h))` via a **shared-denominator** local Padé approximant —
`d` numerator polynomials over **one** common denominator `Q`.  This is
the first place the v0.2 vector pipeline runs end-to-end: it composes
the vector Taylor-jet layer (`VectorCoefficients`) with the shared-`Q`
Hermite–Padé primitive (`SharedPade`).

## Where this fits — the additive v0.2 vector lift

`VectorStepper` is the vector analogue of the scalar `PadeStepper`.  It
mirrors that module's *structure* — jet → `h^k` rescale → Padé →
evaluate at `t = 1` → mutate state — but differs in two essentials:

  1. **First-order, so the state is just `(z, y)`.**  The scalar
     `PadeStepper` integrates a *second*-order ODE `u'' = f(z, u, u')`,
     so it tracks `(z, u, u')` and recovers `u'(z+h)` by analytically
     differentiating the local Padé (chain rule, one extra Horner
     sweep).  A first-order system `y' = f(z, y)` carries no derivative
     in its state: `vector_pade_step!` evaluates only the *value*
     `y(z+h)`, never a Padé derivative.  There is no chain-rule step.

  2. **Shared-`Q` conversion, not scalar `robust_pade`.**  The scalar
     stepper builds one `(m,m)` Padé per (scalar) jet.  The vector
     stepper feeds all `d` rescaled jets to
     `SharedPade.shared_denominator_pade`, which returns `d` numerators
     `P₁,…,P_d` over **one** denominator `Q`.  Component `i`'s new
     value is `Pᵢ(1) / Q(1)` — every component divided by the *same*
     `Q(1)`.  See "Why one shared `Q`" below.

The state object `VectorPadeStepperState{T}` is deliberately minimal —
two fields, `z` and `y::Vector{T}`.  We do not cache the approximant
across calls; each step re-computes from scratch.  Trajectory storage
and dense output are the concern of the higher-level `VectorProblems`
driver (bead V3b), kept out of this layer.

## The 5-step vector algorithm

For `vector_pade_step!(state, f, order, h)`:

  1. **Vector Taylor jet.**  Call `vector_taylor_coefficients(f,
     state.z, state.y, order)` to get `d` length-`(order+1)` coefficient
     vectors — the local Taylor expansion of each component of
     `y(state.z + h')` about `h' = 0`.
  2. **Rescale each jet by `h^k`.**  Form `c̃⁽ⁱ⁾_k = h^k · c⁽ⁱ⁾_k` for
     every component `i` — the Taylor coefficients of `yᵢ(z + h·t)` in
     the variable `t`.  This is the FW 2011 §3.2 variable substitution
     `h' → h·t`; see "Why we rescale by `h^k`" below.  The rescaling is
     applied uniformly to all `d` jets so the shared `Q` lives in the
     same rescaled variable `t`.
  3. **Shared-`Q` Padé.**  Call `shared_denominator_pade(rescaled_jets,
     m)` with the diagonal degree `m = order ÷ 2` (see "Diagonal-Padé
     default" below).  This returns `d` numerator polynomials
     `numerators[i]` and one denominator `denominator`, all in the
     rescaled variable `t` and normalised so `Q(0) = denominator[1] =
     1`.
  4. **Evaluate at `t = 1`.**  Component `i`'s new value is
     `Pᵢ(1) / Q(1)` — Horner on `numerators[i]`, Horner on
     `denominator`, then divide.  All `d` numerators are divided by the
     **same** `Q(1)`.  Because of the `h^k` rescaling, `t = 1`
     corresponds to the original step `h`: the Padé is *local* to the
     expansion centre, in the rescaled variable.
  5. **Mutate state.**  `state.z += h`; `state.y` becomes the new
     `d`-vector.  Return the same `state` object.

Each subcall inherits the fail-fast contract (Rule 1): a too-short jet,
all-zero jets, `Q(0) ≈ 0`, a non-isolated shared null space — all throw
inside `shared_denominator_pade`.  This module adds one further check:
a denominator that is (essentially) zero at `t = 1` — `|Q(1)| ≤ √ε·‖Q‖`,
meaning the step landed on or within roundoff of a pole — throws
`DomainError` with a shorten-the-step suggestion.  The relative `√ε·‖Q‖`
threshold (not bare `iszero`) is needed because a step that lands on a
pole gives `Q(1)` of the order of floating-point roundoff on `‖Q‖`, not
an exact zero; dividing by such a `Q(1)` would yield a meaningless
`~10¹⁵` value rather than a true large solution value.

## Why one shared `Q` — the v0.2 keystone

A vector meromorphic ODE solution's components all blow up at the
**same** movable poles: the singularities are a property of the flow,
not of the individual component (ADR-0019; Mano–Tsuda 2017 §2.2).
Fitting each component with an independent scalar Padé yields `d`
*different* denominators, hence `d` inconsistent pole sets that agree
only up to noise.  The shared `Q` makes the recovered pole set a
single, consistent estimate — the prerequisite for the v0.2 vector
pole-tracking pipeline (bead V7 extracts poles from `denominator`).
Operationally, the shared `Q` also means a single SVD per step rather
than `d` of them.

## Why we rescale by `h^k`

`SharedPade.shared_denominator_pade` (and the scalar `RobustPade` it
generalises) carries scale-sensitive tolerances: the `τ = tol·‖c‖`
singular-value cutoff and the trailing-near-zero trim.  Near a pole the
raw Taylor coefficients of the solution span many orders of magnitude,
and those thresholds misfire on the unscaled jet.  The FW 2011 §3.2 fix
(`references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md` §3.2) absorbs the step length
into the variable: substitute `h' = h·t`, so the series in `h'` becomes
`Σ c_k h^k · t^k = Σ c̃_k t^k` with `c̃_k = h^k · c_k`.  The rescaled
coefficients reflect each term's actual contribution to the value at
`h`, so they have comparable magnitudes for a well-chosen step.  After
Padé-converting `c̃`, evaluation at `t = 1` recovers `Σ c_k h^k`.  This
is the identical rationale documented at length in `PadeStepper.jl`'s
"Why we rescale by `h^k`" section; the vector stepper applies the same
substitution uniformly to all `d` component jets.

## Diagonal-Padé default `m = order ÷ 2`

`shared_denominator_pade` takes a single (shared) denominator degree
`m`.  We mirror the scalar stepper's diagonal choice `m = order ÷ 2`:
for `order = 30` this is the FW 2011 §5.1 `(15,15)` default.  The
shared-`Q` numerators it produces have degree `≤ m` (the upper
Toeplitz block is `(m+1)×(m+1)`), so the approximant is the vector
analogue of a diagonal `(m,m)` Padé.  For odd `order`, integer division
biases `m` downward by one — acceptable, as the FW recipe is
empirically insensitive to ±1 in the degree provided diagonality is
preserved.

## References

  - FW 2011 §3.2 (step rescaling), §5.1 (`(15,15)` default) —
    `references/markdown/FW2011_painleve_methodology_JCP230/
    FW2011_painleve_methodology_JCP230.md`.
  - `src/PadeStepper.jl` — the scalar second-order stepper this module
    mirrors structurally; the `d=1` correctness oracle (VS.1.1).
  - `src/VectorCoefficients.jl` — the vector Taylor-jet layer (ADR-0020).
  - `src/SharedPade.jl` — the shared-`Q` type-II Hermite–Padé primitive
    (ADR-0019); the source of the `m = order÷2` shared denominator.
  - `docs/v0p2_plan.md` — §"The 10 phases", V3a row; §"Architecture
    decision — additive".
  - ADR-0019 (shared-`Q`), ADR-0020 (vector-`y` storage).
"""
module VectorStepper

using LinearAlgebra:        norm
using ..VectorCoefficients: vector_taylor_coefficients
using ..SharedPade:         shared_denominator_pade
using ..SharedPadeDispatch: shared_pade_select

export VectorPadeStepperState, vector_pade_step!, vector_pade_step_with_pade!

# -----------------------------------------------------------------------------
# State
# -----------------------------------------------------------------------------

"""
    VectorPadeStepperState{T}(z, y)

Mutable state for a first-order vector Padé–Taylor stepper integrating
`y' = f(z, y)`, `y ∈ ℂ^d`.  Two fields, element type `T`:

  - `z::T`         — current independent variable.
  - `y::Vector{T}` — `y(z)`, the `d`-component state at the current `z`.

The inner constructor coerces both arguments to `T`:
`VectorPadeStepperState{Float64}(0, [1, 2])` works even though `0`,
`1`, `2` are `Int`s.  Unlike the scalar `PadeStepperState` (three
scalar fields), the state carries no derivative — a first-order system
stores only the value `y`, never `y'`.
"""
mutable struct VectorPadeStepperState{T}
    z::T
    y::Vector{T}

    # Inner constructor: coerce to T so callers can pass mixed-type
    # literals.  `Vector{T}(y)` rebuilds the component vector at element
    # type T; `T(z)` coerces the scalar.  Same recursion-trap avoidance
    # as `PadeStepper.PadeStepperState`.
    VectorPadeStepperState{T}(z, y) where {T} =
        new{T}(T(z), Vector{T}(y))
end

# -----------------------------------------------------------------------------
# One step
# -----------------------------------------------------------------------------

"""
    vector_pade_step!(state::VectorPadeStepperState{T}, f, order::Int,
                      h::Number) -> state

Take one Padé–Taylor step of the first-order vector ODE `y' = f(z, y)`,
mutating `state` in place: `state.z += h`, and `state.y` becomes the
`d`-vector `y(z+h)` produced by the shared-denominator local Padé
approximant in the rescaled variable `t = h'/h`, evaluated at `t = 1`.

`h` is a `Number` rather than a `Real`: when `T <: Complex`, `h` may be
the complex displacement `h·exp(iθ)` along a wedge direction; when
`T <: Real`, `h` must be real (the `T(h)` coercion enforces this at
runtime via `InexactError`).

The 5-step vector algorithm, the shared-`Q` rationale, and the `h^k`
rescaling are documented in the module docstring.

Throws if any subcall throws — the failure modes inherited from
`vector_taylor_coefficients` (`order < 1`, empty `y0`) and
`shared_denominator_pade` (jet too short, all-zero jets, `Q(0) ≈ 0`,
non-isolated null space) — or the local `DomainError` for a shared
denominator that vanishes exactly at `t = 1`.
"""
function vector_pade_step!(state::VectorPadeStepperState{T}, f,
                           order::Int, h::Number) where {T}
    vector_pade_step_with_pade!(state, f, order, h)
    return state
end

"""
    vector_pade_step_with_pade!(state::VectorPadeStepperState{T}, f,
                                order::Int, h::Number)
        -> (state, numerators::Vector{Vector{T}}, denominator::Vector{T})

Like `vector_pade_step!`, but additionally returns the per-step
shared-denominator approximant: the `d` numerator polynomials
`numerators` and the single denominator `denominator`, all in the
*rescaled* variable `t = h'/h` and normalised so `Q(0) = 1`.  The
downstream vector path-network bead (V7) extracts the system's poles
from the roots of `denominator`.

Mirrors `PadeStepper.pade_step_with_pade!`: `vector_pade_step!` is the
wrapper that discards the approximant.  Same fail-fast contract.
"""
function vector_pade_step_with_pade!(state::VectorPadeStepperState{T}, f,
                                     order::Int, h::Number) where {T}
    order ≥ 1 || throw(ArgumentError(
        "vector_pade_step_with_pade!: order must be ≥ 1 (got $order); the " *
        "first-order vector Taylor recursion needs at least one pass."))

    h_T = T(h)

    # Step 1: d Taylor jets of the components about the current z.
    jets = vector_taylor_coefficients(f, state.z, state.y, order)

    # Step 2: rescale every jet, c̃⁽ⁱ⁾_k = h^k · c⁽ⁱ⁾_k, so the shared
    # Padé lives in the unit variable t = h'/h (FW 2011 §3.2).
    rescaled = [_rescale_by_powers(jet, h_T) for jet in jets]

    # Step 3: shared-Q Padé — d numerators over ONE denominator.  The
    # diagonal degree m = order ÷ 2 mirrors the scalar stepper's (m,m).
    # ADR-0028 dispatch: for d ≥ 2, build both the diagonal (m,m) cell A and the
    # wide-square Mano–Tsuda cell B and select per-step by the relative ODE defect
    # (needs f, the current point state.z, and the step h_T — all local here).
    # d = 1 short-circuits to cell A bit-identically (the scalar contract), so the
    # call is uniform.  Returns the same 2-tuple as shared_denominator_pade.
    m = order ÷ 2
    numerators, denominator = shared_pade_select(rescaled, m, f, state.z, h_T)

    # Step 4: evaluate at t = 1.  Q(1) is computed once; every numerator
    # is divided by the SAME Q(1) — the shared-denominator keystone.
    #
    # A step that lands *on* a pole gives Q(1) = 0 in exact arithmetic;
    # in floating point it gives a tiny Q(1) of the order of roundoff on
    # ‖Q‖ (e.g. Q(1) ≈ −2·10⁻¹⁶ for Q = [1, −1] at h = 1).  Dividing by
    # such a Q(1) yields a meaningless ~10¹⁵ "value", not a true large
    # solution value — so we treat |Q(1)| ≤ √ε·‖Q‖ as a pole hit and
    # throw (Rule 1, fail loud), rather than only the exact `iszero`
    # case.  A genuine large-but-finite solution value sits comfortably
    # above this threshold (its Q(1) is O(1), not O(√ε)).
    one_T = one(T)
    q_at_one = _eval_poly(denominator, one_T)
    eps_T = real(T) <: AbstractFloat ? sqrt(eps(real(T))) : sqrt(eps(Float64))
    abs(q_at_one) > eps_T * norm(denominator) || throw(DomainError(h_T,
        "vector_pade_step_with_pade!: shared denominator Q(1) ≈ $q_at_one " *
        "(≤ √ε·‖Q‖) at h=$h; the step landed on (or essentially on) a " *
        "pole of the local shared-Q Padé approximant.  Suggestion: " *
        "shorten the step length, or step around the pole in the " *
        "complex plane."))
    new_y = T[_eval_poly(num, one_T) / q_at_one for num in numerators]

    # Step 5: mutate state.
    state.z = state.z + h_T
    state.y = new_y
    return state, numerators, denominator
end

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

"""
    _rescale_by_powers(c::AbstractVector{T}, h::T) -> Vector{T}

Return the rescaled coefficient vector `c̃[k+1] = h^k · c[k+1]`,
implementing the variable substitution `h' → h·t` (module docstring,
"Why we rescale by `h^k`").  An accumulator `h_pow *= h` avoids a fresh
`h^k` per index.  Identical in spirit to `PadeStepper._rescale_by_powers`;
duplicated here (rather than imported from a sibling module) to keep the
two stepper modules independent of each other's internals.
"""
function _rescale_by_powers(c::AbstractVector{T}, h::T) where {T}
    out = Vector{T}(undef, length(c))
    h_pow = one(T)
    @inbounds for k in 1:length(c)
        out[k] = c[k] * h_pow
        h_pow = h_pow * h
    end
    return out
end

"""
    _eval_poly(c::AbstractVector{T}, t::T) -> T

Evaluate the polynomial `Σ c[k+1] tᵏ` (coefficients low-to-high) by
Horner's rule.  Used to evaluate each shared-`Q` numerator and the
shared denominator at `t = 1`.  A near-zero result on the denominator
is detected and turned into a fail-loud `DomainError` by the caller —
this helper itself does no division and never throws.
"""
function _eval_poly(c::AbstractVector{T}, t::T) where {T}
    s = zero(T)
    @inbounds for k in length(c):-1:1
        s = s * t + c[k]
    end
    return s
end

end # module VectorStepper
