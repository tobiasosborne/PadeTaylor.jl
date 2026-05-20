# ADR-0020 — Vector-y jet storage: `Vector{Taylor1{T}}`

**Status**: Accepted (2026-05-20) | **Bead**: `padetaylor-0ln.7` (V2 —
`VectorCoefficients.jl` vector Taylor jet) | **Plan**: `docs/v0p2_plan.md`
row V2.

## Decision

`VectorCoefficients.vector_taylor_coefficients(f, z0, y0, order)` — the
v0.2 vector Taylor-jet generator for first-order systems `y' = f(z, y)`,
`y ∈ ℂ^d` — carries the *in-flight* vector jet as a **`Vector{Taylor1{T}}`**:
`d` independent scalar `TaylorSeries.jl::Taylor1{T}` series, one per
component. The user-facing **return** value is a `Vector{Vector{T}}` —
`d` plain coefficient vectors `[c_0,…,c_order]`, the per-component shape
the rest of the v0.2 pipeline consumes.

The component-wise bootstrap (FW 2011 §2.1.2 method (b), lifted) advances
all `d` series in lock-step: each pass `f_t = f(z, y)` is a `d`-vector of
`Taylor1`, and `y[i][k] = f_t[i][k-1] / k` for every component `i`. At
`d = 1` this is bit-identical to the scalar `Coefficients.taylor_coefficients_1st`.

## Context

v0.1's `Coefficients.taylor_coefficients_1st` generates the Taylor jet of
a **scalar** ODE solution and stores it as a single `Taylor1{T}` (in
flight) / `Vector{T}` (returned). v0.2 lifts the whole PadeTaylor pipeline
to vector systems — the Noumi–Yamada `A_n⁽¹⁾` higher-order Painlevé
systems and the `P_I⁽²⁾` first-order 4-vector companion form (`docs/
v0p2_pillarB_*` and `docs/v0p2_pillarC_*`). The first design question of
the vector lift is: **how is the `d`-component jet represented?**

Three candidates were considered. The choice propagates through every
later v0.2 module (`VectorStepper`, `VectorProblems`, the vector RHS
factories), so it is recorded here before any of them are built.

## Alternatives considered

### A — `Vector{Taylor1{T}}` (chosen)

`d` independent scalar `Taylor1{T}` series.

- **Composes with `TaylorSeries.jl` for free.** A vector RHS
  `f(z, y) -> ẏ` is just ordinary Julia array + `Taylor1` arithmetic:
  `[y[2], -y[1]]`, `[z^2 + y[1]^2]`, the `P_I⁽²⁾` companion expression
  `−10y[2]^2 − 20y[1]*y[3] − …`. Every operation in the RHS is a scalar
  `Taylor1` op that `TaylorSeries.jl` already implements and that the
  Stage-0 probe (`RESEARCH.md §3.3`) validated for `T ∈ {Float64,
  BigFloat, Arb}`. No custom element type, no operator overloads.
- **Each component hands directly to `SharedPade`.** `y[i].coeffs` is the
  per-component coefficient vector that `SharedPade.shared_denominator_pade`
  takes as one of its `d` input jets (ADR-0019 — the shared-`Q`
  block-Toeplitz stack is built one block per component). The storage
  and the consumer agree with no transposition or repacking.
- **Type-stable.** Each series is seeded `Taylor1(T[y0[i]], order)`;
  `eltype` of every returned vector is exactly `T`, no silent widening
  (the same discipline as the scalar module).
- **`d = 1` reduction is exact.** A length-1 `Vector{Taylor1{T}}` runs
  the identical arithmetic the scalar path runs on its single `Taylor1`,
  so `vector_taylor_coefficients(f, z0, [y0], order)[1]` is bit-identical
  to `taylor_coefficients_1st` — the primary correctness oracle (test
  VC.1.1).

### B — `Taylor1` with vector-valued coefficients (`Taylor1{Vector{T}}`)

A single `Taylor1` whose `k`-th coefficient is the `d`-vector of all
components' `c_k`. **Rejected.** It would require `TaylorSeries.jl` to do
arithmetic over a `Vector{T}` element type — element-wise `*`, `/`, and
the transcendental ops the RHS may use are *not* what `Taylor1{Vector{T}}`
provides (`Taylor1` multiplication is a Cauchy convolution, which is wrong
for a vector-of-components element). The RHS could not be written as
plain array algebra; every system would need a bespoke arithmetic layer —
exactly the hand-rolled coefficient layer CLAUDE.md's hallucination-risk
callout forbids.

### C — a `d × (order+1)` matrix

One dense `Matrix{T}`, row `i` = component `i`'s coefficients.
**Rejected** as the *in-flight* representation: the bootstrap loop and the
RHS still need `Taylor1` arithmetic, so a matrix would be converted to and
from `Vector{Taylor1{T}}` every pass — strictly more work than holding the
`Taylor1` form throughout. The matrix also obscures the per-component
`coeffs` view `SharedPade` wants. (The returned `Vector{Vector{T}}` is
morally a ragged matrix, but as `d` separate vectors it hands to
`SharedPade`'s `jets` argument with no slicing.)

## Consequences

- `VectorCoefficients.jl` is a new module (≤ 200 effective LOC, literate
  top docstring); v0.1 `Coefficients.jl` is untouched and bit-identical —
  the additive architecture of `docs/v0p2_plan.md` §"Architecture
  decision — additive".
- The vector RHS contract is fixed: `f(z, y)` takes a scalar / `Taylor1`
  `z` and an `AbstractVector` `y`, and returns a `d`-vector `ẏ`. The
  later vector RHS factories (Noumi–Yamada, `P_I⁽²⁾`) must produce
  callables of exactly this shape.
- Downstream v0.2 modules (`VectorStepper`, `VectorProblems`) inherit the
  `Vector{Vector{T}}` jet shape and the `Vector{Taylor1{T}}` in-flight
  shape; the shared-`Q` Padé stage (ADR-0019) consumes the returned
  per-component vectors directly.
- Fail-loud (Rule 1): `order < 1` and empty `y0` throw `ArgumentError`
  with informative messages; a RHS that returns a wrong-length vector
  throws inside the bootstrap loop rather than producing a silently
  truncated jet.

## References

- `src/VectorCoefficients.jl` — the implementation (literate, ≤ 200 LOC).
- `test/vector_coefficients_test.jl` — the V2 spec-from-scratch suite
  (158 assertions, GREEN) + the VC.1.6 mutation-proof record.
- `src/Coefficients.jl` — the scalar `taylor_coefficients_1st` this
  module lifts component-wise; the `d = 1` correctness oracle.
- FW 2011 §2.1.2 method (b) —
  `references/markdown/FW2011_painleve_methodology_JCP230/
  FW2011_painleve_methodology_JCP230.md:96-107` — the bootstrap algorithm.
- `external/TaylorIntegration.jl/src/integrator/jetcoeffs.jl` — the
  canonical `Vector{Taylor1}` bootstrap idiom.
- `docs/v0p2_plan.md` — §"Architecture decision — additive", V2 row.
- ADR-0019 — the shared-denominator Padé consumer of these per-component
  jets.
- ADR-0001 — the four-layer architecture this jet layer slots into.
- CLAUDE.md Rules 1, 4, 6, 9, 10.
