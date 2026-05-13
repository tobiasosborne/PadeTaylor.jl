# Worklog 009 — Phase 7: CommonSolveAdapter package extension

**Date**: 2026-05-13
**Author**: Claude Opus
**Scope**: Phase 7 SciML adapter (`PadeTaylorCommonSolveExt`) per
ADR-0003.  Translation layer wiring `(PadeTaylorProblem,
PadeTaylorAlg)` into `CommonSolve.jl`'s `init` / `step!` / `solve!` /
`solve` interface.  374/374 tests GREEN (340 prior + 34 new across 7
testsets).  Bead `padetaylor-2vz` closes.

> Take-home: a deliberately minimal extension — `PadeTaylorAlg` lives
> in the main module so users get the un-qualified name after `using
> CommonSolve`; the extension only adds the four `CommonSolve` methods.
> ADR-0003's "Translation only — no algorithmic logic in the extension"
> contract is the literal scope.  Trajectory bytes are bit-identical to
> `solve_pade`.

## Public-API shape

```julia
using PadeTaylor, CommonSolve

prob = PadeTaylorProblem(f, (u0, up0), zspan; order = 30)
alg  = PadeTaylorAlg(; h_max = 0.5, max_steps = 100_000)

# Convenience: CommonSolve's default `solve = solve! ∘ init`.
sol = solve(prob, alg)

# Streaming form:
integ = init(prob, alg)
while !integ.done
    step!(integ)
end
sol = solve!(integ)
```

The bit-identity contract: `solve(prob, alg)` returns a
`PadeTaylorSolution` whose `z, y, h, pade` vectors are `==` to
`solve_pade(prob; h_max = alg.h_max, max_steps = alg.max_steps)`.
Tested by CS.1.1 / CS.1.3.

## Design — why `PadeTaylorAlg` lives in the main module

ADR-0003 specifies `PadeTaylorAlg <: CommonSolve.AbstractAlgorithm`
inside the extension.  In practice this forces users to write the
qualified name `PadeTaylor.PadeTaylorCommonSolveExt.PadeTaylorAlg(...)`
to construct an algorithm — clunky.

**v1 amendment**: move the struct declaration to `src/PadeTaylor.jl`
and drop the `AbstractAlgorithm` subtyping.  `CommonSolve.jl` uses
**duck-typing** (`solve(prob, alg)` dispatches on the concrete pair,
not on `<: AbstractAlgorithm`), so this is functionally equivalent
and gives ergonomic access:

```julia
using PadeTaylor, CommonSolve
alg = PadeTaylorAlg(; h_max = 0.5)    # no qualified path needed
```

The extension still owns the methods `init`, `step!`, `solve!` — those
require `CommonSolve` loaded, which is the actual gating condition.
The struct itself is a plain value type with no extension-dependent
fields, so it costs nothing to declare unconditionally in main.

## What changed

`Project.toml`: added `CommonSolve` to `[extras]` + `[targets].test`
so the test environment resolves it.  (The `[weakdeps]` +
`[extensions]` entries were already present from Phase Z.)

`src/PadeTaylor.jl`: declared `PadeTaylorAlg{H<:Real}` struct + kwarg
constructor + `export`.  ~15 LOC including docstring.

`ext/PadeTaylorCommonSolveExt.jl`: 130 LOC including the literate
docstring.  Defines:
  - `PadeTaylorIntegrator{F, T, Y, P, H}` — mutable streaming state.
  - `CommonSolve.init(prob, alg)` — constructs the integrator,
    pre-pushes IC.  Validates `h_max > 0` and 2nd-order branch.
  - `CommonSolve.step!(integ)` — one `pade_step_with_pade!` per call,
    with `h_step = min(alg.h_max, z_end - state.z)` so the final step
    lands EXACTLY on `zspan[2]` (Mutation G target).  No-op on done
    integrator (Mutation H target).
  - `CommonSolve.solve!(integ)` — loop step! until done; return
    `PadeTaylorSolution` from accumulators.

`test/ext_commonsolve_test.jl`: 7 testsets + mutation-proof commentary.
  - CS.1.1: `solve ≡ solve_pade` (bit-identical trajectory).
  - CS.1.2: streaming `init + step!` loop ≡ direct solve.
  - CS.1.3: `solve!` after manual init ≡ solve.
  - CS.2.1: streaming — uneven final step, done flag, no-op-on-done.
  - CS.2.2: degenerate zspan handled at construction layer.
  - CS.3.1: fail-fast on h_max ≤ 0.
  - CS.3.2: fail-fast on max_steps overrun.

## Mutation-proof (verified 2026-05-13)

**Mutation G** — in `step!`, drop the `min(h_max, z_end - state.z)`
clamp: `h_step = h_max_T`.  This breaks the "final step lands exactly
on `z_end`" invariant.

Verified bite: 2 fails on CS.2.1 — `integ.z_vec[end] ≈ 0.55` evaluates
`0.6 ≈ 0.55` (overshoot); `integ.h_vec[end] ≈ 0.15` evaluates `0.2 ≈
0.15` (unclamped step).  Initial probe of this mutation **did not
bite** because CS.1.1 used h_max=1.5 on zspan=(0, 1.5) — a single
even step.  The fix was to upgrade CS.2.1 to use zspan=(0, 0.55) with
h_max=0.2, forcing the final 0.15 step to differ from h_max.  Lesson:
**tests with evenly-divisible final steps cannot bite the clamping logic.**

**Mutation H** — in `step!`, drop the `integ.done && return integ`
early-out.  This breaks the "step! after done is a no-op" contract.

Verified bite: 1 fail on CS.2.1 — `length(integ.z_vec)` after the
extra `step!` on a done integrator is 5 instead of 4.  The inner
`pade_step_with_pade!` accepts `h_step = 0` and silently appends a
zero-length segment, so the mutation doesn't throw — it just corrupts
the accumulator.  Caught by the no-op assertion.

Both mutations restored before commit per CLAUDE.md Rule 4.

## Frictions surfaced

  - **F1. `Pkg.resolve()` fails on existing manifest**: the repo's
    `Manifest.toml` has `GenericLinearAlgebra v0.4` but `Project.toml`
    pins `compat = "0.3"`.  `Pkg.test()` re-resolves internally with
    downgrade to 0.3.19 (works), but `Pkg.resolve()` standalone errors
    on the version pin conflict.  Not blocking; existing project state.
  - **F2. Initial Project.toml edit didn't appear to land**: my first
    `Edit` of `[extras]` reported success but `cat`'ing the file
    showed unchanged content.  Re-applied; second edit took.  Possibly
    a transient file-state cache issue.  Recorded for future suspicion.
  - **F3. Naïve mutation didn't bite**: as noted under Mutation G,
    the chosen test setup happened to have an even final step.  Had to
    redesign CS.2.1 to use uneven `zspan`/`h_max` before mutation G's
    bite was visible.  Same pattern as worklog 008 §"F2 PN.3.1's
    original near-IC framing was trivial" — tests need to exercise
    the code path their mutations target.

## Bead state

Closed in this session:
  - `padetaylor-2vz` — Phase 7 CommonSolveAdapter (GREEN at this commit).

Still open from prior sessions (no change):
  - `padetaylor-rgp`, `padetaylor-c2p`, `padetaylor-k31`,
    `padetaylor-kvi`, `padetaylor-jhq`, `padetaylor-bvh`,
    `padetaylor-grc`, `padetaylor-61j`, `padetaylor-8pi`.

## Pointers

  - [`Project.toml`](../../Project.toml) — `[weakdeps]` + `[extensions]`.
  - [`src/PadeTaylor.jl`](../../src/PadeTaylor.jl) — `PadeTaylorAlg`
    declaration.
  - [`ext/PadeTaylorCommonSolveExt.jl`](../../ext/PadeTaylorCommonSolveExt.jl)
    — the extension methods.
  - [`test/ext_commonsolve_test.jl`](../../test/ext_commonsolve_test.jl)
    — 7 testsets + mutation-proof commentary.
  - [`docs/adr/0003-extensions-pattern.md`](../adr/0003-extensions-pattern.md)
    — extensions design (with v1 amendment about `PadeTaylorAlg`
    location).
  - [CommonSolve.jl](https://github.com/SciML/CommonSolve.jl) source —
    25-LOC interface; we add four methods.
