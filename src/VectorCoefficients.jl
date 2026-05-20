"""
    PadeTaylor.VectorCoefficients

Generate the local Taylor jet of a first-order **vector** ODE solution
`y' = f(z, y)`, `y ∈ ℂ^d`, at an expansion centre `z0`, returning a
`Vector{Vector{T}}` — `d` coefficient vectors, one per component, each
of the form `[c_0, c_1, …, c_order]`.

## What this module is for

v0.1's `Coefficients.taylor_coefficients_1st` produces the Taylor jet of
a **scalar** ODE solution.  v0.2 lifts the whole PadeTaylor pipeline from
scalar analytic ODEs to vector systems — the Noumi–Yamada higher-order
Painlevé systems `A_n⁽¹⁾`, and the first Painlevé hierarchy member
`P_I⁽²⁾` (a first-order 4-vector companion system).  This module is the
vector analogue of the v0.1 jet layer: the "give me a high-order Taylor
jet of the *vector* solution" stage that everything downstream consumes.

It is deliberately **RHS-agnostic**.  It takes any callable
`f(z, y) -> ẏ` where `z` is a scalar / `Taylor1` and `y` is a `d`-vector,
and returns the jet of every component.  The concrete vector RHS
factories (Noumi–Yamada `A_n⁽¹⁾`, the `P_I⁽²⁾` 4-vector form) are built
on top of it in later v0.2 beads — this module knows nothing about them.

## The additive-architecture rationale

The v0.2 vector lift is **additive** (v0.2 plan, §"Architecture decision
— additive"): new modules alongside the v0.1 scalar modules, which stay
bit-identical.  `VectorCoefficients` is a *new* module; `Coefficients.jl`
is untouched.  This keeps the v0.1 scalar tests bit-identical (backward
compatibility is non-negotiable) and respects the Rule 6 ≤200-LOC cap by
not bolting a vector code path onto the scalar file.

## The algorithm — FW 2011 §2.1.2 method (b), lifted component-wise

The scalar bootstrap (FW 2011 §2.1.2 method (b);
`references/markdown/FW2011_painleve_methodology_JCP230/
FW2011_painleve_methodology_JCP230.md:96-107`) loops on `Taylor1{T}`
arithmetic: seed `y` with its constant term, then for `k = 1,…,order`
evaluate `f_t = f(z, y)` and set `y[k] = f_t[k-1] / k`.  Each pass gains
exactly one new correct coefficient.

For a vector ODE this lifts *component-wise* with no change of idea.  The
unknown `y` becomes a `d`-vector of `Taylor1{T}` — one independent series
per component.  Each pass `f_t = f(z, y)` returns a `d`-vector of
`Taylor1{T}` (the RHS, applied to the partially-known vector jet).  The
integration relation `y' = f` holds component-wise, so

    y[i][k] = f_t[i][k-1] / k        for every component i = 1,…,d.

The correctness argument is unchanged: on entry to pass `k`, every
component's first `k` coefficients are correct, hence each `f_t[i][k-1]`
is correct, hence each `y[i][k]` set this pass is correct.  This is
exactly the shape of `TaylorIntegration.jl::jetcoeffs!`
(`external/TaylorIntegration.jl/src/integrator/jetcoeffs.jl`), which
advances a `Vector{Taylor1}` one order at a time inside its integrator;
here we expose it as a standalone jet generator.

## The d = 1 reduction — the primary correctness oracle

For `d = 1` the loop sets `y[1][k] = f_t[1][k-1] / k`, the *same*
arithmetic the scalar `taylor_coefficients_1st` performs on its single
`Taylor1`.  Hence `vector_taylor_coefficients(f, z0, [y0], order)[1]` is
**bit-identical** to `taylor_coefficients_1st(scalar_f, z0, y0, order)`
for the corresponding scalar RHS.  Test VC.1.1 asserts exactly this; it
is the cheapest, tightest correctness check available.

## The `Vector{Taylor1{T}}` representation

The vector jet is carried as a `Vector{Taylor1{T}}` — `d` independent
`Taylor1` series.  ADR-0020 records why this beats the alternatives
(a `Taylor1` with vector-valued coefficients; a `d×(order+1)` matrix):
it composes with `TaylorSeries.jl`'s scalar `Taylor1` arithmetic
component-wise for free (`f(z, y)` is just ordinary array+Taylor algebra,
no custom element type), and each component's `coeffs` vector hands
directly to `SharedPade.shared_denominator_pade`'s per-component jet
input.

## Type-stability watch

`y0::AbstractVector{T}` is *not* silently widened.  Each component's seed
is `Taylor1(T[y0[i]], order)` — an explicitly `T`-typed coefficient
vector — and `z = z0 + Taylor1(T, order)` keeps the same `T` flowing
through the loop.  The element type of every returned vector is exactly
`T`.  Test VC.1.4 asserts this for `T = BigFloat` (256-bit).

## References

  - FW 2011 §2.1.2 method (b) —
    `references/markdown/FW2011_painleve_methodology_JCP230/
    FW2011_painleve_methodology_JCP230.md:96-107`.
  - `src/Coefficients.jl` — the scalar `taylor_coefficients_1st` this
    module lifts component-wise; the v0.1 jet layer and the d=1 oracle.
  - `external/TaylorIntegration.jl/src/integrator/jetcoeffs.jl` — the
    canonical `Vector{Taylor1}` bootstrap idiom.
  - `docs/v0p2_plan.md` — §"Architecture decision — additive", V2 row.
  - `docs/adr/0020-vector-y-storage.md` — the `Vector{Taylor1{T}}`
    storage decision.
  - ADR-0001 (own-repo) — the four-layer architecture this slots into.
"""
module VectorCoefficients

using TaylorSeries: Taylor1

export vector_taylor_coefficients

# -----------------------------------------------------------------------------
# First-order vector ODE: y' = f(z, y), y ∈ ℂ^d
# -----------------------------------------------------------------------------

"""
    vector_taylor_coefficients(f, z0::T, y0::AbstractVector{T},
                               order::Integer) where {T}
        -> Vector{Vector{T}}

For the first-order vector ODE `y' = f(z, y)` with initial condition
`y(z0) = y0` (a `d`-vector), return the Taylor jet of the solution as a
`Vector{Vector{T}}`: `d` coefficient vectors, where the `i`-th vector is
`[c_0, c_1, …, c_order]` for component `i` of `y(z0 + h)` about `h = 0`.

`f(z, y)` must take a scalar / `Taylor1` `z` and a `d`-vector `y` and
return the `d`-vector `ẏ`.  Each returned vector has length `order + 1`
and element type exactly `T` (no silent widening).

Bootstrap-by-substitution per FW 2011 §2.1.2 method (b), lifted
component-wise; see the module docstring for the algorithm, the `d = 1`
reduction to `Coefficients.taylor_coefficients_1st`, and citations.

Throws `ArgumentError` (Rule 1 — fail loud) for:

  - `order < 1` — a zero-order "expansion" is just `y0` and carries no
    information about `f`; we surface this as an input error rather than
    silently returning `[y0]`-vectors.
  - empty `y0` — a zero-component system has nothing to integrate.
"""
function vector_taylor_coefficients(f, z0::T, y0::AbstractVector{T},
                                    order::Integer) where {T}
    order < 1 && throw(ArgumentError(
        "vector_taylor_coefficients requires order ≥ 1 (got $order); a " *
        "zero-order expansion carries no information about f."))

    d = length(y0)
    d == 0 && throw(ArgumentError(
        "vector_taylor_coefficients requires a non-empty y0; got an " *
        "empty initial-condition vector (d = 0 — nothing to integrate)."))

    # Symbolic expansion variable: z = z0 + h, with h = Taylor1(T, order).
    z = z0 + Taylor1(T, order)

    # Seed each component's series with its constant term y0[i] only.
    # Explicit `T[y0[i]]` makes the underlying coefficient vector
    # Vector{T}, not Vector{Float64} — mirrors the scalar module.
    y = Taylor1{T}[Taylor1(T[y0[i]], order) for i in 1:d]

    # Bootstrap loop: each pass gains exactly one new correct coefficient
    # in every component.  On entry to pass k, every component's first k
    # coefficients are correct, so each f_t[i][k-1] is correct, and the
    # integration relation y' = f gives y[i][k] = f_t[i][k-1] / k.
    @inbounds for k in 1:order
        f_t = f(z, y)
        length(f_t) == d || throw(ArgumentError(
            "vector_taylor_coefficients: f(z, y) returned a length-" *
            "$(length(f_t)) vector but y0 has d = $d components; the " *
            "RHS must map a d-vector to a d-vector."))
        for i in 1:d
            y[i][k] = f_t[i][k-1] / k
        end
    end

    # `Taylor1.coeffs` is a `FixedSizeVectorDefault{T}` (TaylorSeries
    # v0.20+); the public API contract is `Vector{T}`, so copy out.
    return Vector{T}[Vector{T}(y[i].coeffs) for i in 1:d]
end

end # module VectorCoefficients
