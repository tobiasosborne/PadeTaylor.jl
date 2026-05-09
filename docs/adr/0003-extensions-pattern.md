# ADR-0003 — Extensions pattern for Arb and CommonSolve

**Status**: Accepted (2026-05-09)
**Context**: Stage 1 design lock. The PRD requires `Arblib.jl`-backed
arb-prec as first-class for v1. The design also calls for an optional
`CommonSolve.jl` adapter so PadeTaylor.jl can participate in SciML's
`solve(prob, alg; ...)` ecosystem. Neither should be a hard dep —
they balloon the load time and, in `Arblib`'s case, drag in Flint's
C-FFI surface.

## Decision

Use Julia 1.9+ **package extensions** (a.k.a. weak deps) for both:

```toml
[deps]
TaylorSeries          = "..."
GenericLinearAlgebra  = "..."
Polynomials           = "..."
LinearAlgebra         = "..."   # stdlib
Printf                = "..."   # stdlib

[weakdeps]
Arblib                = "..."
CommonSolve           = "..."

[extensions]
PadeTaylorArblibExt    = "Arblib"
PadeTaylorCommonSolveExt = "CommonSolve"
```

Extension files at `ext/PadeTaylorArblibExt.jl` and
`ext/PadeTaylorCommonSolveExt.jl`. Each is a single-file Julia
module that imports the parent package's internals via
`PadeTaylor: pade_svd, default_tol, ...` and adds methods.

## Why extensions

- **Core dep graph minimal.** `using PadeTaylor` brings in
  `TaylorSeries`, `GenericLinearAlgebra`, `Polynomials` — three
  small, pure-Julia packages. No FFI, no Flint, no SciML stack.
- **Opt-in features.** A user who only needs `Float64` integration
  doesn't pay for `Arblib`'s precompile time. A user who doesn't
  use SciML doesn't pay for `CommonSolve`'s API surface.
- **Canonical Julia idiom (1.9+).** Packages like
  `OrdinaryDiffEq.jl`, `ChainRulesCore.jl`, and most of SciML use
  this pattern.
- **Load triggered automatically.** `using PadeTaylor; using Arblib`
  loads the extension; the user does not need to import anything
  from `PadeTaylor.LinAlg` directly.

## What goes in each extension

### `PadeTaylorArblibExt.jl` (~80 LOC)

Methods to add when `Arblib` is loaded:
- `LinAlg.pade_svd(A::Matrix{Arblib.Arb})` — convert to
  `Matrix{BigFloat}` (mid-point rounding, radius discarded; see
  ADR-0002 caveat), dispatch to `GenericLinearAlgebra.svd`.
- `LinAlg.pade_svd(A::Matrix{Arblib.Acb})` — same, via `BigComplex`.
- `default_tol(::Type{Arblib.Arb})` — `Arb`-precision-aware tolerance
  `2.0^(-precision(BigFloat) + 10)`.
- `Coefficients.taylor_coefficients_*` extra paths if `Arb` exposes
  features `BigFloat` does not (e.g. ball radius reporting). Prefer
  default `T <: Number` paths first; specialise only where measurably
  needed.

### `PadeTaylorCommonSolveExt.jl` (~50 LOC)

- `struct PadeTaylorAlg <: CommonSolve.AbstractAlgorithm; opts; end`
- `CommonSolve.solve(prob::PadeTaylorProblem, ::PadeTaylorAlg;
  kwargs...) = solve_pade(prob; merge(alg.opts, kwargs)...)`
- `CommonSolve.init`, `step!`, `solve!` for streaming use.

Translation only — no algorithmic logic in the extension.

## Versioning constraint

Julia 1.9+ is required for the extensions pattern. Document this in
`Project.toml`'s `[compat]` section: `julia = "1.9"`.

We adopt **Julia 1.10** as the actual minimum supported version
(LTS-track; matches scientist-workbench's "Bun 1.3+" baseline
discipline of choosing a stable, current release).

## Citations

- `RESEARCH.md §7.1 Q5` — API decision for SciML adapter.
- `RESEARCH.md §5.1` and ADR-0002 — `Arblib` route through
  `BigFloat`.
- Julia docs on package extensions (any 1.9+ release notes).
- `OrdinaryDiffEq.jl` and `ChainRulesCore.jl` extensions as reference
  patterns.

## Consequences

- `Project.toml` has the `[weakdeps]` and `[extensions]` sections
  per Julia 1.9 spec.
- Two ext files at `ext/`. Tests for them go under
  `test/ext/<name>_test.jl` and run conditionally on the extension
  being loaded.
- Documentation will note: "for arb-prec, also `using Arblib`; for
  SciML integration, also `using CommonSolve`".

## Alternatives considered, rejected

- **Hard deps on `Arblib` + `CommonSolve`.** Rejected: balloons load
  time; doesn't reflect actual user needs.
- **Re-export `pade_svd` from a separate `PadeTaylorArb` package.**
  Rejected: forces users to know about a second package; fragments
  the API.
- **Conditional `if isdefined(Main, :Arblib)` checks.** Rejected:
  brittle, runtime-dispatch, doesn't compose with precompilation.
