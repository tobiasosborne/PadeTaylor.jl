# Worklog 003 — Phase 5 (PadeStepper) and the inner-loop integration

**Date**: 2026-05-09
**Author**: Claude (orchestrating); opus subagent (implementation)
**Scope**: Phase 5 of Stage 2 — `PadeStepper`. The big integration:
the four-layer architecture's first end-to-end run.

> Take-home from this shard: **trust the test suite, listen to the
> subagent's empirical findings**.  Two algorithmic choices that
> looked obvious from the spec turned out wrong on the actual numerical
> test cases; the subagent surfaced both with concrete error numbers,
> and the docstring-as-chapter records the empirical justification.

## What changed

`src/PadeStepper.jl` (~87 LOC body + 145-line literate docstring):

- `mutable struct PadeStepperState{T}` with fields `z, u, up`.
- `pade_step!(state, f, order, h) -> state` — one step of the 2nd-order
  Padé-Taylor pipeline.
- Private helpers: `_rescale_by_powers`, `_evaluate_pade`,
  `_evaluate_pade_deriv` (quotient rule).

`test/padestepper_test.jl`: 4 testsets / 16 `@test` calls. Pinned
oracle in `test/_oracle_padestepper.jl` (three-source: Mathematica
closed-form ℘ ≡ Mathematica NDSolve ≡ mpmath.odefun at 40 dps).

`external/probes/padestepper-oracle/` — capture.{wl, py} + verify.jl.

`DESIGN.md §4 Phase 5` — corrected to match the shipped 5-step
algorithm + acknowledge the three test deviations from the original
table.

## Phase-5 acceptance: 207/207 GREEN

## Two algorithmic deviations from the original DESIGN.md sketch

### Deviation 1 — `h^k` rescaling of Taylor coefs before robust_pade

DESIGN's sketch went straight from `coefs_u = taylor_coefficients_2nd(...)`
to `P_u = robust_pade(coefs_u, m, n)`. That fails catastrophically
on test 5.1.4.

The reason is concrete and worth recording: `RobustPade.robust_pade`
faithfully ports Chebfun's `padeapprox.m`, which carries two scale-
sensitive thresholds. Specifically the "function ≡ 0" early-exit
fires when `max|c[1:m+1]| ≤ tol · max|c|`. At `z = 0.9` near a
Weierstrass-℘ pole at `z = 1`, the raw Taylor coefs of `u(z + h')`
about `h' = 0` have `c_0 ≈ 100` and `c_30 ≈ 3·10³³` — coefficients
that grow as `1/(distance to pole)^k` per the Cauchy estimate. With
`tol = 1e-14`, `tol · ‖c‖∞ ≈ 3·10¹⁹`, comfortably above `c_0 = 100`.
The early-exit fires; `robust_pade` returns the literal zero
function; `pade_step!` returns garbage.

The standard fix (FW 2011 §3.2 line 396) is to absorb the step length
into the variable: `h' = h · t`, so the Taylor series becomes
`Σ c_k h^k · t^k = Σ c̃_k t^k`. The rescaled coefficients `c̃_k`
reflect the *actual contribution to the value at `h`*, so they have
comparable magnitudes whenever the step is well-chosen, and the
threshold no longer misfires. After Padé-converting `c̃`, evaluation
at `t = 1` recovers the original `Σ c_k h^k`.

The subagent surfaced this empirically: tried the spec's no-rescale
approach first, hit the 5.1.4 failure, traced through `robust_pade`'s
internals, found the threshold misfire, added the rescaling.
Recorded in the module docstring as "Why we rescale by h^k".

This is the kind of finding that confirms the spec discipline: the
DESIGN sketch was correct *in principle* (per the FW algorithm), but
the implementation detail of *where* the rescaling lives (in
PadeStepper, not in RobustPade) only crystallised once we saw
empirically what failed.

### Deviation 2 — analytic differentiation of P_u beats re-Padé for u'

DESIGN's sketch had Phase 5 do *two* `robust_pade` calls: one for the
Taylor coefs of `u`, one for the formally-differentiated coefs of
`u'`. The structural argument: the analytic derivative of an `(m, n)`
rational is `(m + n - 1, 2n)`, breaking diagonal-Padé balance.

Empirically, the re-Padé path loses ~40× accuracy on test 5.1.2
(compose two ℘-steps near the pole): re-Padé gives `up` rel-err
`2.7·10⁻¹¹` vs `7·10⁻¹³` for analytic differentiation. The
explanation: the formally-differentiated coefficients `(k+1)·c_{k+1}`
have a *steeper* growth than `c_k` itself by the factor `k+1`, so
after `h^k` rescaling the resulting `c̃'` is more skewed than `c̃`.
The `(14, 14)` Padé from `c̃'` inherits worse conditioning than the
parent `(15, 15)` from `c̃`. The structural diagonal-Padé argument
*is* relevant to extrapolation across a *range* of points; we
evaluate at exactly one point (`t = 1`), so the degree-imbalance
argument does not apply.

The subagent's choice: skip the second `robust_pade` call, evaluate
`P_u'(1) / h` via the quotient rule on the *parent* Padé. Saves one
SVD per step (~30% wall time) and 40× accuracy improvement on the
near-pole composition test. Recorded in the docstring as "Why
analytic differentiation beats re-Padé for u'".

## Frictions surfaced

### F1. Inner-vs-outer constructor recursion trap

The subagent's first attempt:
```julia
mutable struct PadeStepperState{T}
    z::T; u::T; up::T
end
PadeStepperState{T}(z, u, up) where {T} =
    PadeStepperState{T}(T(z), T(u), T(up))
```

This `StackOverflowError`s. The default inner constructor
`PadeStepperState{T}(z::T, u::T, up::T)` requires *exactly typed* args.
The outer constructor `PadeStepperState{T}(z, u, up)` matches anything,
including the `T(z), T(u), T(up)` call inside it — so it dispatches
to itself, infinite loop.

Fix: replace with an inner constructor:
```julia
mutable struct PadeStepperState{T}
    z::T; u::T; up::T
    PadeStepperState{T}(z, u, up) where {T} = new{T}(T(z), T(u), T(up))
end
```

The `new{T}(...)` form bypasses dispatch, so it always reaches the
type-bound construction. **Lesson**: Julia inner constructors are
the right tool when you want to coerce-then-allocate; the outer-
constructor-with-recursive-dispatch trick used in many tutorials only
works when the inner constructor is *more restrictive* than the
outer.

### F2. wolframscript variable names with underscores parse as patterns

When writing `capture.wl`, names like `u_05` and `up_05` were parsed
by Mathematica as `Pattern[u, Blank[]]` (followed by a stray `_05`).
Output got nonsense. Fix: use camelCase Mathematica-side
(`uHalf`, `upHalfFromNd`, etc.); keep underscored names only in the
emitted Julia output strings.

### F3. mpmath has no native Weierstrass-P

For ℘ values, the Python cross-check has to rely on `mpmath.odefun`
(arbitrary-precision Taylor-method ODE solver), integrating the same
ODE `u'' = 6u^2` numerically. Three-source consensus on ℘ values is
therefore (closed-form Mathematica) ≡ (Mathematica NDSolve) ≡
(mpmath.odefun); for PI we have only two (no closed form).

### F4. DESIGN.md's `z = 1.36` pole was incorrect

DESIGN test 5.1.4 mentions a "pole at z ≈ 1.363". On the lattice
specified by FW (`c₁ = -1, c₂ = g₃ = 2`), the pole nearest z = 0 is
at `z = -c₁ = 1` (where `WeierstrassP`'s argument is 0).
Wolframscript confirmed: `1/u(1) = 0` exactly. We retargeted 5.1.4
at `z = 0.9 → 0.95` (still 0.05 from the pole, exercising the
near-pole behaviour DESIGN intended).

### F5. DESIGN.md's tritronquée IC for 5.1.3 had no paper-pinned value

We don't have a tritronquée-IC paper locally with a 16-digit pin;
substituted the FW ICs and asserted PI's RHS difference from ℘ as a
sanity check on the `+z` contribution. Test pinned via NDSolve at
WorkingPrecision=50.

## Mutation-proof procedure (verified before commit)

  - **Mutation A** — sign flip on Padé eval point for u
    (`one_T → -one_T` in line 224 of `src/PadeStepper.jl`):
    fails the `state.u` assertion in all 4 testsets.
  - **Mutation C** — flip chain-rule h division on the u' path
    (`/ h_T → * h_T` in line 225):
    fails the `state.up` assertion in all 4 testsets.

Both mutations bite both real-axis-clean (5.1.1, 5.1.3) and
near-pole (5.1.2, 5.1.4) cases, so the test corpus genuinely
exercises the implementation in both regimes.

## Pointers

- [`src/PadeStepper.jl`](../../src/PadeStepper.jl) — module + the
  empirical-rationale docstring.
- [`test/padestepper_test.jl`](../../test/padestepper_test.jl) — 4
  testsets + verified mutation-proof block.
- [`external/probes/padestepper-oracle/`](../../external/probes/padestepper-oracle/)
  — capture (Mathematica closed-form + NDSolve, mpmath.odefun) +
  verifier.
- [`docs/worklog/002-phase-4-spec-correction.md`](002-phase-4-spec-correction.md)
  — previous shard, for context on the broader spec-discipline pattern.
