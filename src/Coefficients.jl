"""
    PadeTaylor.Coefficients

Generate the local Taylor expansion of an ODE solution at a chosen
expansion centre `z0`, returning a plain `Vector{T}` of coefficients
`[c_0, c_1, …, c_order]` ready for hand-off to `RobustPade.robust_pade`
in the next stage of the FW 2011 pipeline.

## What this module is for

This is the *first* of the four-layer architecture (ADR-0001): the
"give me a high-order Taylor jet of the solution" layer.  Everything
downstream (`RobustPade`, `StepControl`, `PadeStepper`) consumes its
output as a coefficient vector and treats the ODE as opaque — only
this module ever calls back into the user-supplied right-hand side.

The two functions exported here, `taylor_coefficients_1st` and
`taylor_coefficients_2nd`, cover first- and second-order scalar ODEs;
that is the full surface area we need for the FW 2011 reproduction
target (Painlevé I, equianharmonic Weierstrass-℘, etc.).  Higher-order
or vector ODEs are deliberately out of scope for v1.

## The algorithm — FW 2011 §2.1.2 method (b)

The classical "differentiate the ODE n times" approach (FW 2011 §2.1.1,
`references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md:88-94`) explodes combinatorially
and becomes unworkable beyond order ~10.  Method (b)
(`...md:96-107`) trades that explosion for a clean fixed-point loop
on Taylor1 arithmetic:

    "...at first, use only one term of (2.3), substitute into (2.4)
     and integrate, giving two correct terms.  For each similar step,
     one gains one correct coefficient in the expansion."
                                                      — FW 2011 §2.1.2

Concretely for `dy/dz = f(z, y)`:

  1. Build `z_taylor = z0 + Taylor1(T, order)` — the symbolic `z0 + h`.
  2. Seed `y_taylor = Taylor1(T[y0], order)` — constant term only.
  3. For `k = 1, 2, …, order`:
       evaluate `f_t = f(z_taylor, y_taylor)`.  Because `y_taylor`'s
       first `k` coefficients are already correct, `f_t[k-1]` is also
       correct.  The integration relation `y' = f` ⇒ `y_k = f_{k-1}/k`
       gives one new correct coefficient: `y_taylor[k] = f_t[k-1]/k`.
  4. Return `y_taylor.coeffs`.

For `u'' = f(z, u, u')` we evolve `u` and its formal derivative `up`
in lock-step.  The relation `up_k = (k+1) · u_{k+1}` keeps `up` synced
with `u`; the integration relation `u'' = f` ⇒
`u_j = f_{j-2} / (j · (j-1))` advances `u` two coefficients at a time.

## Why we delegate Taylor arithmetic to `TaylorSeries.jl`

The Stage-0 empirical probe at `external/probes/taylorseries-arb/`
established (see `RESEARCH.md §3.3`) that `TaylorSeries.jl::Taylor1{T}`
correctly handles all the arithmetic and transcendental operations we
need for `T ∈ {Float64, BigFloat, Arb}` at orders up to 80, including
256-bit precision propagation with measured radii below 2⁻²⁵⁸.  Hand-
rolling a coefficient layer would duplicate code, miss optimisations
(`TaylorSeries.jl` uses Horner-style accumulators internally), and
sacrifice the empirical validation already in hand.  CLAUDE.md's
hallucination-risk callout is explicit on this point: do not re-derive
generic Taylor arithmetic.

## Type-stability watch

`y0::T` is *not* silently widened.  We construct the seed `Taylor1` as
`Taylor1(T[y0], order)` (an explicitly `T`-typed coefficient vector),
not `Taylor1([y0], order)` — the latter can default to `Vector{Float64}`
if the literal `y0` is parsed as `Float64`.  The same `T` flows through
`z0 + Taylor1(T, order)` and the whole loop, so the returned vector
satisfies `eltype(out) == T`.  Tests 3.1.1 and 3.1.4 assert this for
`T = Float64` and `T = BigFloat` respectively.

## References

  - FW 2011 §2.1.2 method (b) —
    `references/markdown/FW2011_painleve_methodology_JCP230/
    FW2011_painleve_methodology_JCP230.md:96-107`.
  - `RESEARCH.md §3.3` — `TaylorSeries.jl::Taylor1{T}` empirical
    validation across `T ∈ {Float64, BigFloat, Arb}`.
  - `external/probes/taylorseries-arb/probe.jl` — concrete usage.
  - `external/TaylorIntegration.jl/src/integrator/jetcoeffs.jl:23-43`
    — canonical Julia pattern; we adapt the same shape but expose it
    as a standalone function rather than a step inside an integrator.
  - ADR-0001 (own-repo) — four-layer architecture rationale.
"""
module Coefficients

using TaylorSeries: Taylor1

export taylor_coefficients_1st, taylor_coefficients_2nd

# -----------------------------------------------------------------------------
# First-order: dy/dz = f(z, y)
# -----------------------------------------------------------------------------

"""
    taylor_coefficients_1st(f, z0::T, y0::T, order::Int) where {T} -> Vector{T}

For the first-order ODE `dy/dz = f(z, y)` with initial condition
`y(z0) = y0`, return the Taylor coefficients `[c_0, c_1, …, c_order]`
of `y(z0 + h)` expanded about `h = 0`.  The returned vector has length
`order + 1` and element type exactly `T` (no silent widening).

Bootstrap-by-substitution per FW 2011 §2.1.2 method (b); see the module
docstring for the algorithm and citations.

Throws `ArgumentError` for `order < 1` (a zero-order "expansion" is
just `[y0]` and offers no information about `f`; we surface this as an
input error rather than silently returning `[y0]`).
"""
function taylor_coefficients_1st(f, z0::T, y0::T, order::Int) where {T}
    order < 1 && throw(ArgumentError(
        "taylor_coefficients_1st requires order ≥ 1 (got $order); a " *
        "zero-order expansion carries no information about f."))

    # Symbolic expansion variable: z = z0 + h, with h = Taylor1(T, order).
    z = z0 + Taylor1(T, order)

    # Seed y with the constant term y0 only.  Explicit T[y0] ensures
    # the underlying coefficient vector is Vector{T}, not Vector{Float64}.
    y = Taylor1(T[y0], order)

    # Bootstrap loop: each pass gains exactly one new correct coefficient.
    @inbounds for k in 1:order
        f_t = f(z, y)
        y[k] = f_t[k-1] / k
    end

    # `Taylor1.coeffs` is a `FixedSizeVectorDefault{T}` (TaylorSeries v0.20+);
    # the public API contract is `Vector{T}`, so we copy out.
    return Vector{T}(y.coeffs)
end

# -----------------------------------------------------------------------------
# Second-order: u'' = f(z, u, u')
# -----------------------------------------------------------------------------

"""
    taylor_coefficients_2nd(f, z0::T, y0::T, y1::T, order::Int) where {T} -> Vector{T}

For the second-order ODE `u'' = f(z, u, u')` with initial conditions
`u(z0) = y0`, `u'(z0) = y1`, return the Taylor coefficients
`[c_0, c_1, …, c_order]` of `u(z0 + h)` about `h = 0`.  The returned
vector has length `order + 1` and element type exactly `T`.

Algorithm: evolve `u` and its formal derivative `up = u'` jointly,
maintaining the invariant on entry to pass `j` (j = 2, 3, …, order):

  - `u[0], …, u[j-1]` are correct;
  - `up[0], …, up[j-2]` are correct, related to `u` by the formal-
    differentiation identity `up[k] = (k+1) · u[k+1]`.

At pass `j`:

  - evaluate `f_t = f(z, u, up)`; the coefficient `f_t[j-2]` depends
    only on `u[0..j-2]` and `up[0..j-2]`, both of which are correct;
  - apply the integration relation `u'' = f` ⇒
    `u[j] = f_t[j-2] / (j · (j-1))` — one new correct `u` coefficient;
  - resync `up[j-1] = j · u[j]` to restore the invariant for pass j+1.

The j=2 entry needs no `up` resync at the top of the loop because
`up[0] = y1` is the seeded initial condition.

Throws `ArgumentError` for `order < 2`: with `order = 0` the result is
just `[y0]` (no `f`), with `order = 1` it is `[y0, y1]` (still no `f`).
We surface either as an input error rather than returning silently.
"""
function taylor_coefficients_2nd(f, z0::T, y0::T, y1::T, order::Int) where {T}
    order < 2 && throw(ArgumentError(
        "taylor_coefficients_2nd requires order ≥ 2 (got $order); a " *
        "lower-order expansion carries no information about f."))

    # Symbolic expansion variable.
    z = z0 + Taylor1(T, order)

    # Seed u with [y0, y1, 0, 0, …] (length order+1) and up with the
    # corresponding [y1, 0, 0, …]; both explicitly T-typed so we don't
    # widen.  up's index-0 coefficient is u'(z0) = y1; higher indices
    # are filled in lock-step with u below.
    u  = Taylor1(T[y0, y1], order)
    up = Taylor1(T[y1],     order)

    # Bootstrap loop: each pass gains one new correct u-coefficient.
    # Invariant on entry to pass j (j ≥ 2):
    #   - u_0, …, u_{j-1} are correct;
    #   - up_0, …, up_{j-2} are correct (synced via up_k = (k+1)·u_{k+1}).
    # That is exactly what f_t[j-2] depends on, so the read below is
    # correct, and the integration relation u'' = f then gives one new
    # u coefficient: u_j = f_{j-2} / (j·(j-1)).
    # After setting u[j] we resync up at index j-1 to maintain the
    # invariant for pass j+1.  The j=2 entry is a no-op for up because
    # up_0 = y1 is the seeded initial condition.
    @inbounds for j in 2:order
        f_t = f(z, u, up)
        u[j] = f_t[j-2] / (j * (j-1))
        up[j-1] = j * u[j]
    end

    return Vector{T}(u.coeffs)
end

end # module Coefficients
