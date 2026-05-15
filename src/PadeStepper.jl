"""
    PadeTaylor.PadeStepper

Take one Padé-Taylor step of a second-order analytic ODE
`u'' = f(z, u, u')` from a current state `(z, u, u')` to a new state
`(z+h, u(z+h), u'(z+h))` via the Fornberg–Weideman 2011 inner-loop
recipe: build a high-order local Taylor jet, rescale-and-convert it to
a robust Padé rational, evaluate the rational (and its derivative) at
the corresponding unit point.  This module is the first place in the
four-layer architecture where the whole pipeline runs end-to-end.

## Where this fits — the four-layer architecture (ADR-0001)

The Phase-3 `Coefficients` module produces a coefficient vector
`[c_0, c_1, …, c_order]` of `u(z + h')`'s local Taylor expansion about
the current `z`.  The Phase-2 `RobustPade` module turns that vector
into a `PadeApproximant{T}` — a numerator/denominator pair whose
analytic continuation domain extends *past* the radius of convergence
of the Taylor series, capturing nearby poles of the solution by
arranging for the denominator to vanish at the right place rather than
allowing the numerator to blow up.  This module composes the two: it
is the consumer of both, and the producer for the higher-level
`Problems.solve_pade` driver in Phase 6.

The state object `PadeStepperState{T}` is deliberately minimal — three
scalar fields `z`, `u`, `up`.  We do not cache the `PadeApproximant`
across calls; each step re-computes from scratch.  Storage and dense
output (the *trajectory*) is Phase 6's concern, kept out of this layer.

## The 5-step algorithm

For `pade_step!(state, f, order, h)`:

  1. **Taylor jet.**  Call `Coefficients.taylor_coefficients_2nd(f,
     state.z, state.u, state.up, order)` to get `coefs_u`, the length-
     `(order+1)` Taylor coefficients of `u(z + h')` about `h' = 0`.
  2. **Rescale & Padé of u.**  Form `c̃_k = h^k · c_k` (the Taylor
     coefficients of `u(z + h·t)` in the variable `t`), then build the
     diagonal Padé approximant `P_u = robust_pade(c̃, m, n)` with
     `m = n = order ÷ 2`.  See "Why we rescale by h^k" below; the
     rescaling is essential to keep `RobustPade`'s special-case-zero
     and trim-near-zero thresholds from misfiring on near-pole jets.
  3. **Evaluate at `t = 1`.**  `new_u = P_u(1)` — Horner-Horner-divide.
     Because of the `h^k` rescaling, `t = 1` corresponds to the
     original step `h`, not the original `z + h`: the Padé is *local*
     to the expansion centre, in the rescaled variable.
  4. **Evaluate `u'` via the chain rule.**  `u(z + h') = P_u(h'/h)`
     ⇒ `u'(z + h) = P_u'(1) / h`, where `P_u'` is the analytic
     derivative of the rational `P_u` (numerator-and-denominator
     differentiation; see `_evaluate_pade_deriv`).  No second SVD;
     no second `robust_pade` call.  See "Why analytic differentiation
     beats re-Padé" below for the empirical numbers that motivated
     this choice.
  5. **Mutate state.**  `state.z += h`; `state.u`, `state.up` get the
     new values.  Return the same `state` object.

Each subcall (Coefficients, RobustPade, evaluation) inherits the fail-
fast contract: a singular `C̃` in the SVD step, a denominator that
vanishes at `t = 1`, an unsupported `order` — all throw with context
rather than returning silently-wrong values.

## Why we rescale by `h^k`

`RobustPade.robust_pade` carries two scale-sensitive thresholds (it
faithfully ports Chebfun's `padeapprox.m`, lines 69 and 134):

  - The "function ≡ 0" special case fires when
    `max|c[1:m+1]| ≤ tol · max|c|`.
  - The trailing-near-zero numerator trim fires when
    `|a[k]| ≤ tol · ‖c‖₂`.

Near a pole, the raw Taylor coefficients of the *solution* span
many orders of magnitude — `c_0 ≈ u(z)` may be `O(10²)` while `c_30 ≈
3·10³³` (case 5.1.4 in `test/padestepper_test.jl`, with `z = 0.9` near
a Weierstrass-℘ pole at `z = 1`).  `tol · ‖c‖∞` is then `tol ·
3·10³³ ≈ 3·10¹⁹`, comfortably above `c_0 = 100` — the special-case-
zero check fires spuriously and `robust_pade` returns the literal
function `r ≡ 0`.

The standard FW-2011 fix (FW 2011 §3.2 line 396 in
`references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md`) is to absorb the step length
into the variable: substitute `h' = h · t`, so the Taylor series in
`h'` becomes `Σ c_k h^k · t^k = Σ c̃_k t^k`, with `c̃_k = h^k · c_k`.
The rescaled `c̃` reflect *the coefficients' actual contributions to
the value at `h`*, so they have comparable magnitudes whenever the
step is well-chosen, and the special-case threshold no longer
misfires.  After Padé-converting `c̃`, evaluation at `t = 1` recovers
the original `Σ c_k h^k`.

This is also why the RobustPade module *itself* does not scale: that
layer is a faithful port of `padeapprox.m` (its determinism contract
is bit-identical to the Octave oracle).  The scaling is the
caller's job and lives here.

## Why analytic differentiation beats re-Padé for `u'`

A natural-looking alternative is to build a *separate* Padé from the
formally-differentiated coefficients `coefs_up[k] = (k+1)·coefs_u[k+1]`
and evaluate that.  The argument for it is structural: the analytic
derivative of an `(m, n)` rational is an `(m + n - 1, 2n)` rational,
breaking diagonal Padé's degree-balance, while a fresh Padé of the
derivative coefficients restores diagonality.

Empirically (see the four oracle cases in `test/padestepper_test.jl`),
that argument loses to a more practical effect: the formally-
differentiated coefficients `(k+1)·c_{k+1}` have a *steeper* growth
than `c_k` itself by the factor `k+1`, so the resulting `c̃'` after
`h^k` rescaling is more skewed than `c̃`, and the `(14, 14)` Padé
inherits worse conditioning than the `(15, 15)` parent.  In test
5.1.2 (compose two ℘-steps to land near the pole) the re-Padé path
hits relative error `2.7·10⁻¹¹` for `u'`, while the analytic-
derivative path hits `7·10⁻¹³` — a `40×` improvement, and the
difference between failing and passing the spec.

The analytic-derivative path also saves one `robust_pade` call
(one fewer SVD per step, ≈ 30% wall time on Float64 at order 30).
The argument we lose — that the `(m+n-1, 2n)` rational degrades
extrapolation near the pole boundary — is real, but only matters
for *evaluating P_u' at points distant from `t = 1`*.  We evaluate
at exactly one point (`t = 1`), so the degree-imbalance argument
does not apply.

## Diagonal-Padé default `(order ÷ 2, order ÷ 2)`

FW 2011 §5.1 line 277 in `references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md:271-281` settles on `(15, 15)`
for `order = 30`.  We generalise to `(order÷2, order÷2)`.  For odd
`order`, integer division biases `m = n` downward by one — losing one
denominator-degree; this is acceptable because the FW recipe is
empirically insensitive to ±1 in `(m, n)` provided diagonality is
preserved (RESEARCH.md §7.1).

## Adaptive Padé step size — FFW 2017 §2.1.2

FFW 2017 §2.1.2 (`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
FFW2017_painleve_riemann_surfaces_preprint.md:74-97`) introduces an
adaptive step-size controller specialised to Padé-Taylor stepping.
The estimator is the one-extra-Taylor-coefficient truncation-error
formula `T(h) = |ε_{n+1} h^{n+1} / a(h)|`, where `a(h)` is the
**numerator** polynomial of the local Padé evaluated at the step,
`b_r` are the denominator coefficients (with `b_0 = 1`), and
`ε_k = c_k + Σ_{r=1..ν} b_r c_{k-r}` are the Padé error coefficients
of the partial-sum residual `q(h) w(ζ+h) - p(h) = Σ ε_k h^k` (md:76).

The controller's rescale law (md:88-91 eq. 2) is

    q = (k · Tol / T(h))^(1/(n+1)),    h := q · h,

with the FFW-recommended conservative factor `k = 1e-3`.  Acceptance
is `T(h) ≤ Tol`; on accept, the *next* step's initial `h` is seeded
from the accepted `|q·h|` (md:93) — the controller has memory.

Three helpers ship here:

  - `ffw_truncation_error(f, z, u, up, n, h)` — compute `T(h)` from
    one extra Taylor coefficient.  Internally builds the order-`(n+1)`
    Taylor jet (one pass beyond the `n` used by the step) so that
    `c_{n+1}` is available; the Padé denominator `b` and numerator
    `a` are reused from the rescaled order-`n` jet.  All algebra is
    in the *rescaled* variable `t = h'/h`: `c̃_k = h^k c_k`,
    `ε̃_{n+1} = c̃_{n+1} + Σ b_r c̃_{n+1-r}`, `a(h) ≡ a_rescaled(t=1)`.
    Returns the real-typed magnitude `|ε̃_{n+1} / a_rescaled(1)|`,
    which equals FFW's `|ε_{n+1} h^{n+1} / a(h)|` by construction.

  - `ffw_rescale_q(Tol, T_h, n; k = 1e-3)` — the rescale factor.
    Pure: takes the controller magnitudes and returns `q`.  Safe at
    `T_h = 0` (returns `Inf`, callers must check `T_h > Tol` before
    rescaling).

  - `adaptive_pade_step!(state, f, n, h_init; adaptive_tol, ...)` —
    the full controller: takes the wedge-direction `h_init` (complex
    or real), iterates `T(h) > Tol ⇒ h := q·h` up to `max_rescales`
    times, then runs the accepted step via `pade_step_with_pade!`.
    Returns `(state, P_u, meta)` where `meta` is a NamedTuple of
    `(:h_used, :h_step, :T_h, :n_rescales)`.  Throws
    `ErrorException` if `max_rescales` is exceeded (Rule 1).

The controller is *complex-h aware*: `h_init` carries the wedge-step
direction (e.g. `h_step = h·exp(i·θ)` from PathNetwork's wedge), and
rescaling shrinks the magnitude only — multiplying a complex `h` by
a real `q ∈ (0, 1]` preserves direction.  `meta.h_used` reports the
real magnitude of the accepted step (used by `PathNetwork` to seed
the next step's initial `h`).

## References

  - FW 2011 §3.1 (pole handling rationale), §3.2 (step rescaling),
    §5.1 line 277 (`(15, 15)` default) —
    `references/markdown/FW2011_painleve_methodology_JCP230/
    FW2011_painleve_methodology_JCP230.md:271-281`.
  - FFW 2017 §2.1.2 (adaptive Padé step) —
    `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/
    FFW2017_painleve_riemann_surfaces_preprint.md:74-97`.
  - `src/Coefficients.jl` — Phase-3 Taylor jet generator.
  - `src/RobustPade.jl` — Phase-2 `PadeApproximant{T}` + `robust_pade`.
  - `external/probes/padestepper-oracle/` — three-source oracle
    (Mathematica `WeierstrassP` + `NDSolve`, `mpmath.odefun`).
  - ADR-0001 — four-layer architecture rationale.
  - ADR-0011 — adaptive Padé step.
"""
module PadeStepper

using ..Coefficients: taylor_coefficients_2nd
using ..RobustPade:   PadeApproximant, robust_pade

export PadeStepperState, pade_step!, pade_step_with_pade!
export ffw_truncation_error, ffw_rescale_q, adaptive_pade_step!

# -----------------------------------------------------------------------------
# State
# -----------------------------------------------------------------------------

"""
    PadeStepperState{T}(z, u, up)

Mutable state for one Padé-Taylor stepper integrating a 2nd-order ODE
`u'' = f(z, u, u')`.  Three scalar fields, all of element type `T`:

  - `z::T`  — current independent variable.
  - `u::T`  — `u(z)` at the current point.
  - `up::T` — `u'(z)` at the current point.

The constructor coerces its arguments to `T` (`PadeStepperState{Float64}(0,
1, 2)` works even though `0`, `1`, `2` are `Int`s).
"""
mutable struct PadeStepperState{T}
    z::T
    u::T
    up::T

    # Inner constructor: coerce to T so callers can pass mixed-type
    # literals (`PadeStepperState{Float64}(0, 1.5, 2)` works) without
    # the outer-constructor / default-inner-constructor recursion trap
    # (an outer `PadeStepperState{T}(z, u, up) where {T}` form would
    # call itself, since computed `T(z)` values are dispatch-equivalent
    # to the unconstrained outer signature).
    PadeStepperState{T}(z, u, up) where {T} = new{T}(T(z), T(u), T(up))
end

# -----------------------------------------------------------------------------
# One step
# -----------------------------------------------------------------------------

"""
    pade_step!(state::PadeStepperState{T}, f, order::Int, h::Number) -> state

Take one Padé-Taylor step of the 2nd-order ODE `u'' = f(z, u, u')`,
mutating `state` in place: `state.z += h`, and `state.u`, `state.up`
become the values of `u(z+h)`, `u'(z+h)` produced by the local diagonal
Padé approximation in the rescaled variable `t = h'/h`.

`h` is a `Number` rather than a `Real`: when `T <: Complex` (path-
network mode), `h` is the complex displacement `h·exp(iθ)` along a
wedge direction.  When `T <: Real`, `h` must be real — the inner
`T(h)` coercion enforces this at runtime via `InexactError`.

The 5-step algorithm and the rationale for `h^k` rescaling, diagonal
Padé, evaluation at `t = 1`, and analytic differentiation for `u'`
are in the module docstring.

Throws if any subcall throws — the failure modes inherited from
`Coefficients` (`order < 2`) and `RobustPade` (singular SVD, near-zero
denominator), or the local check for a denominator that vanishes
exactly at `t = 1` (`DomainError`).
"""
function pade_step!(state::PadeStepperState{T}, f,
                    order::Int, h::Number) where {T}
    pade_step_with_pade!(state, f, order, h)
    return state
end

"""
    pade_step_with_pade!(state::PadeStepperState{T}, f, order::Int, h::Number)
        -> (state, P_u::PadeApproximant{T})

Like `pade_step!`, but additionally returns the per-step Padé
approximant `P_u` built from the rescaled Taylor coefficients
`c̃_k = h^k · c_k`.  `P_u` is the rational approximant in the
*rescaled* variable `t = h'/h`: `u(state.z_old + h'·) = P_u(h'/h)`,
so it interpolates the segment `[state.z_old, state.z_old + h]`
when evaluated at `t ∈ [0, 1]`.

This is the function `Problems.solve_pade` calls inside its main
loop: it needs both the advanced state (to keep stepping) *and*
the Padé (to store for dense interpolation, per the per-segment
Padé idea in `Problems.PadeTaylorSolution`).  `pade_step!` is the
backwards-compatible wrapper that discards the Padé — Phase 5's
existing tests do not need the Padé and continue to work
unchanged.

`PathNetwork.path_network_solve` (Phase 10) also calls this, passing
a *complex* `h` representing a wedge-direction step in the complex
plane.  In that mode, `T = Complex{S}` and the rescaled coefficients
`c̃_k = h^k · c_k` are complex; `RobustPade.robust_pade` handles the
complex linear algebra unchanged.

Same fail-fast contract as `pade_step!`.
"""
function pade_step_with_pade!(state::PadeStepperState{T}, f,
                              order::Int, h::Number) where {T}
    order ≥ 2 || throw(ArgumentError(
        "pade_step_with_pade!: order must be ≥ 2 (got $order); the " *
        "2nd-order Taylor recursion needs at least two passes."))

    h_T = T(h)

    # Step 1: Taylor coefficients of u about the current z.
    coefs_u = taylor_coefficients_2nd(f, state.z, state.u, state.up, order)

    # Step 2: rescale c̃_k = h^k · c_k and Padé in the unit variable.
    coefs_u_scaled = _rescale_by_powers(coefs_u, h_T)
    m = n = order ÷ 2
    P_u = robust_pade(coefs_u_scaled, m, n)

    # Steps 3 & 4: evaluate u and u' at t = 1.  u'(z+h) = P_u'(1)/h
    # follows from the chain rule on u(z + h') = P_u(h'/h).
    one_T  = one(T)
    new_u  = _evaluate_pade(P_u, one_T)
    new_up = _evaluate_pade_deriv(P_u, one_T) / h_T

    # Step 5: mutate state.
    state.z  = state.z + h_T
    state.u  = new_u
    state.up = new_up
    return state, P_u
end

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

"""
    _rescale_by_powers(c::Vector{T}, h::T) -> Vector{T}

Return the rescaled coefficient vector `c̃[k+1] = h^k · c[k+1]`,
implementing the variable substitution `h' → h · t` discussed in the
module docstring.  Computed by an accumulator `h_pow *= h` to avoid
the `O(N)` cost of a fresh `h^k` per coefficient and the catastrophic
loss of precision at large `k` from a separate `pow` per index.
"""
function _rescale_by_powers(c::Vector{T}, h::T) where {T}
    out = Vector{T}(undef, length(c))
    h_pow = one(T)
    @inbounds for k in 1:length(c)
        out[k] = c[k] * h_pow
        h_pow = h_pow * h
    end
    return out
end

"""
    _evaluate_pade(P::PadeApproximant{T}, z::T) -> T

Evaluate `P(z) = (Σ a[k+1] z^k) / (Σ b[k+1] z^k)` with Horner's rule
on numerator and denominator separately, then divide.  Throws
`DomainError(z, …)` if the denominator vanishes exactly at `z` —
treated as a fail-fast signal that the step landed on a pole of the
local Padé (Rule 1).

Two-pass Horner (rather than fused) keeps the numerical paths
identical to a hand evaluation of `Σ a[k+1] z^k`; we don't try to
balance numerator-denominator scale here because `b` is already
normalised to `b[1] = 1` by `RobustPade._trim_and_normalise`, so the
denominator's leading magnitude is `O(1)` for small `z`.
"""
function _evaluate_pade(P::PadeApproximant{T}, z::T) where {T}
    num = zero(T)
    @inbounds for k in length(P.a):-1:1
        num = num * z + P.a[k]
    end
    den = zero(T)
    @inbounds for k in length(P.b):-1:1
        den = den * z + P.b[k]
    end
    iszero(den) && throw(DomainError(z,
        "PadeStepper._evaluate_pade: denominator vanishes at z=$z; the " *
        "step landed exactly on a pole of the local Padé approximant. " *
        "Suggestion: shorten the step length, or step around the pole " *
        "in the complex plane."))
    return num / den
end

"""
    _evaluate_pade_deriv(P::PadeApproximant{T}, z::T) -> T

Evaluate `(d/dz) P(z)` via the quotient rule:
`P' = (N' · D - N · D') / D²`, where `N(z) = Σ a[k+1] z^k`,
`D(z) = Σ b[k+1] z^k`, and the derivative coefficients are
`N'[k] = k · a[k+1]` (and similarly for `D'`).

Implemented by four Horner sweeps (N, N', D, D') sharing the loop
counter; cost is one `_evaluate_pade` plus three Horner passes.
Throws `DomainError(z, …)` if `D(z) = 0` (same fail-fast pole-
detection as `_evaluate_pade`).
"""
function _evaluate_pade_deriv(P::PadeApproximant{T}, z::T) where {T}
    n = length(P.a)
    m = length(P.b)
    N = zero(T); @inbounds for k in n:-1:1; N = N * z + P.a[k]; end
    D = zero(T); @inbounds for k in m:-1:1; D = D * z + P.b[k]; end
    Nt = zero(T); @inbounds for k in n:-1:2; Nt = Nt * z + (k - 1) * P.a[k]; end
    Dt = zero(T); @inbounds for k in m:-1:2; Dt = Dt * z + (k - 1) * P.b[k]; end
    iszero(D) && throw(DomainError(z,
        "PadeStepper._evaluate_pade_deriv: denominator vanishes at z=$z; " *
        "the step landed exactly on a pole of the local Padé approximant. " *
        "Suggestion: shorten the step length."))
    return (Nt * D - N * Dt) / (D * D)
end

# -----------------------------------------------------------------------------
# Adaptive Padé step — FFW 2017 §2.1.2
# -----------------------------------------------------------------------------

"""
    ffw_truncation_error(f, z::T, u::T, up::T, order::Int, h::Number)
        -> Real

Compute the FFW 2017 §2.1.2 truncation-error estimate
`T(h) = |ε_{n+1} · h^{n+1} / a(h)|` for a single tentative Padé step
of size `h` from the current state `(z, u, up)` of the second-order
ODE `u'' = f(z, u, up)`.  `order = n` is the Taylor truncation order
(typically 30); the local Padé is `(ν, ν)` with `ν = n ÷ 2`.

Internally we build the order-`(n+1)` Taylor jet so that the extra
coefficient `c_{n+1}` is available (FFW md:74).  The Padé is built
from the rescaled order-`n` slice `c̃_k = h^k · c_k` per the FW 2011
rescaling discipline (see module docstring "Why we rescale by h^k").
The Padé error coefficient is then evaluated in the *rescaled*
variable:

    ε̃_{n+1} = c̃_{n+1} + Σ_{r=1..ν} b_r · c̃_{n+1-r},

and the controller magnitude is `|ε̃_{n+1} / a_rescaled(1)|`, which
is algebraically identical to FFW's `|ε_{n+1} h^{n+1} / a(h)|` after
substituting `c̃_k = h^k c_k` and `a_rescaled(1) = a(h)`.

`h` may be complex (wedge-direction step); we return a *real-typed*
magnitude.  The denominator `a_rescaled(1) = 0` case throws
`DomainError` per Rule 1: a Padé numerator that vanishes at the step
endpoint is a pathological singularity, not a small `T(h)` value.

Throws `ArgumentError` if `order < 2` (the 2nd-order recursion needs
at least two passes, and the `n+1` Taylor expansion needs `order+1 ≥
3` so the `ε_{n+1}` formula has a valid lookup target).
"""
function ffw_truncation_error(f, z::T, u::T, up::T,
                              order::Int, h::Number) where {T}
    order ≥ 2 || throw(ArgumentError(
        "ffw_truncation_error: order must be ≥ 2 (got $order)."))

    h_T = T(h)

    # Build order-(n+1) Taylor so c_{n+1} is in coefs_u[end].
    coefs_u = taylor_coefficients_2nd(f, z, u, up, order + 1)

    # Rescale c̃_k = h^k · c_k for k = 0..n+1 (length n+2).
    coefs_u_scaled = _rescale_by_powers(coefs_u, h_T)

    # Build the diagonal (ν, ν) Padé from the rescaled order-n slice
    # (first n+1 entries, matching pade_step_with_pade!'s convention).
    n_int = order
    m = n_int ÷ 2
    coefs_n  = @view coefs_u_scaled[1:(n_int + 1)]
    P_u = robust_pade(Vector{T}(coefs_n), m, m)

    # ε̃_{n+1} = c̃_{n+1} + Σ_{r=1..ν} b_r · c̃_{n+1-r}.  P_u.b is
    # `[b_0, b_1, …, b_ν_eff]` with `b_0 = 1`; `ν_eff` may be < ν if
    # `_trim_and_normalise` dropped trailing near-zero denominator
    # entries.  The Σ runs up to `r = ν_eff` only — terms with
    # `r > ν_eff` correspond to `b_r = 0` and contribute nothing.
    cn_plus_1 = coefs_u_scaled[n_int + 2]    # c̃_{n+1} (1-based: index n+2).
    eps_nplus1 = cn_plus_1
    νeff = length(P_u.b) - 1
    @inbounds for r in 1:νeff
        eps_nplus1 += P_u.b[r + 1] * coefs_u_scaled[n_int + 2 - r]
    end

    # Numerator a(h) evaluated at t = 1 in the rescaled variable.
    a_at_one = _eval_poly_at_one(P_u.a)
    iszero(a_at_one) && throw(DomainError(h,
        "ffw_truncation_error: rescaled Padé numerator a_rescaled(1) = 0 " *
        "(equivalently a(h) = 0) at h=$h; the step landed on a zero of " *
        "the local Padé.  Suggestion: shorten h or shift step direction."))

    return abs(eps_nplus1 / a_at_one)
end

"""
    ffw_rescale_q(Tol::Real, T_h::Real, order::Int; k::Real = 1.0e-3)
        -> Real

Compute the FFW 2017 §2.1.2 rescale factor

    q = (k · Tol / T_h)^(1/(n+1))

per eq. (2) at md:88-91.  `k = 1e-3` is FFW's recommended conservative
factor.  Returns `Inf` when `T_h = 0` (the step is already exact at
working precision; the caller's `T_h ≤ Tol` accept-test fires before
`q` is computed in normal use, so this branch is only reached if the
caller asks for `q` at `T_h = 0` directly — defensive, no throw).

Pure; no state, no allocations beyond the return scalar.
"""
function ffw_rescale_q(Tol::Real, T_h::Real, order::Int; k::Real = 1.0e-3)
    Tol > 0 || throw(ArgumentError(
        "ffw_rescale_q: Tol must be positive (got $Tol)."))
    T_h ≥ 0 || throw(ArgumentError(
        "ffw_rescale_q: T_h must be non-negative (got $T_h)."))
    k > 0 || throw(ArgumentError(
        "ffw_rescale_q: k must be positive (got $k)."))
    order ≥ 1 || throw(ArgumentError(
        "ffw_rescale_q: order must be ≥ 1 (got $order)."))
    T_h == 0 && return Inf
    return (k * Tol / T_h) ^ (1 / (order + 1))
end

"""
    adaptive_pade_step!(state::PadeStepperState{T}, f, order::Int,
                        h_init::Number;
                        adaptive_tol::Real = 1.0e-12,
                        k_conservative::Real = 1.0e-3,
                        max_rescales::Int = 50)
        -> (state, P_u, meta::NamedTuple)

Take one **adaptive** Padé-Taylor step.  Iterate `T(h) > Tol ⇒ h := q·h`
(FFW 2017 §2.1.2 eq. 2 rescale) until either the truncation-error
estimate satisfies `T(h) ≤ Tol` or `max_rescales` rescales have been
attempted (Rule 1 fail-loud — does NOT silently accept a too-large
step).  On accept, run the standard `pade_step_with_pade!` with the
accepted `h` and return the new state, its Padé, and a metadata
NamedTuple.

`h_init` may be complex (wedge-direction step); rescaling multiplies
by a real `q ∈ (0, 1]`, preserving the wedge direction.

`meta`:
  - `h_used::Real`     — magnitude of the accepted step `|h|`.
  - `h_step::T`        — the *complex* accepted step `h_init·∏q_k`.
  - `T_h::Real`        — final truncation-error estimate (≤ Tol).
  - `n_rescales::Int`  — number of rescale iterations consumed.

The caller (typically `PathNetwork.path_network_solve`) uses
`meta.h_used` as the *seed* `h` for the next step (FFW md:93 "the
initial step length is always the scaled step length stored at the
current point").

Throws:
  - `ErrorException` if `max_rescales` exhausted with `T(h) > Tol`.
  - `ArgumentError` if `adaptive_tol ≤ 0`, `k_conservative ≤ 0`, or
    `max_rescales < 1`.
"""
function adaptive_pade_step!(state::PadeStepperState{T}, f, order::Int,
                             h_init::Number;
                             adaptive_tol::Real = 1.0e-12,
                             k_conservative::Real = 1.0e-3,
                             max_rescales::Int = 50) where {T}
    adaptive_tol > 0 || throw(ArgumentError(
        "adaptive_pade_step!: adaptive_tol must be positive (got $adaptive_tol)."))
    k_conservative > 0 || throw(ArgumentError(
        "adaptive_pade_step!: k_conservative must be positive (got $k_conservative)."))
    max_rescales ≥ 1 || throw(ArgumentError(
        "adaptive_pade_step!: max_rescales must be ≥ 1 (got $max_rescales)."))

    h = T(h_init)
    n_rescales = 0
    Th = ffw_truncation_error(f, state.z, state.u, state.up, order, h)

    while Th > adaptive_tol
        n_rescales ≥ max_rescales && throw(ErrorException(
            "adaptive_pade_step!: max_rescales=$max_rescales exhausted at " *
            "z=$(state.z) with T(h)=$Th > Tol=$adaptive_tol.  Suggestion: " *
            "tighten the controller (smaller h_init, larger max_rescales) " *
            "or relax adaptive_tol."))
        q = ffw_rescale_q(adaptive_tol, Th, order; k = k_conservative)
        h *= q
        n_rescales += 1
        Th = ffw_truncation_error(f, state.z, state.u, state.up, order, h)
    end

    # Accepted: run the actual step.
    _, P_u = pade_step_with_pade!(state, f, order, h)
    meta = (h_used     = abs(h),
            h_step     = h,
            T_h        = Th,
            n_rescales = n_rescales)
    return state, P_u, meta
end

# Evaluate the polynomial `a` (low-to-high) at t = 1 — a Horner-style
# sum.  Used by `ffw_truncation_error` to compute `a_rescaled(1)`.
# Kept private; not exported.
function _eval_poly_at_one(a::AbstractVector{T}) where {T}
    s = zero(T)
    @inbounds for k in eachindex(a)
        s += a[k]
    end
    return s
end

end # module PadeStepper
