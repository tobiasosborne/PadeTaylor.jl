# Worklog 005 — Phase 6 (Problems) implementation: the Padé-vs-Taylor pole-bridge demo

**Date**: 2026-05-09
**Author**: Claude (orchestrating); opus subagent (implementation)
**Scope**: Phase 6 of Stage 2 — `Problems` module, shipped per the
revised acceptance criteria of `docs/worklog/004-phase-6-pivot.md`.
This shard is the GREEN report for the pivoted scope: `solve_pade` +
`PadeTaylorProblem` + `PadeTaylorSolution` plus a `taylor_eval` helper
for the side-by-side Padé-vs-Taylor headline test.  It also records
one algorithmic finding caught at impl time (the order/rtol coupling
in test 6.1.6), in the lineage of worklog 002's spec-drift pattern.

> Take-home from this shard: at `z = 1.05` (just past the lattice
> pole at `z = 1`), the same-Taylor-jet Padé conversion beats plain
> truncation by **≈ 9.86 orders of magnitude** on the Fornberg-
> Weideman 2011 test problem — the analytic-continuation advantage of
> the algorithm is now empirically demonstrable in one test setup.
> The cost: tightening the BF-256 test (6.1.6) revealed that the
> `(15, 15)` Padé approximation-error floor at `t = 0.7` dominates
> Float64-precision noise at order 30; only at order 40 does
> BF-256 add value over Float64.

## Orchestration pattern

Phase 6 used a serial three-step pattern over a single Julia process
(per CLAUDE.md Rule 7 — no parallel Julia agents):

  1. **Orchestrator writes RED tests + extends the oracle.**  Edited
     `external/probes/problems-oracle/{capture.wl, capture.py,
     verify.jl}` to pin the new pole-bridge points; replaced the
     placeholder `test/problems_test.jl` with the 6 testsets from
     `docs/worklog/004-phase-6-pivot.md`'s test plan.
  2. **Subagent implements GREEN.**  Sole `julia` invocation owner
     for the run; produces `src/Problems.jl` and updates the umbrella
     re-export list in `src/PadeTaylor.jl`.
  3. **Orchestrator verifies + mutation-proves.**  Re-runs `Pkg.test()`
     to confirm 218/218 GREEN, applies both planned mutations,
     observes RED, restores.

Read-only literature subagents (FW 2011, GGT 2013 line-range checks)
ran in parallel where useful; nothing held a `julia` lock.

## What changed

### 1. `external/probes/problems-oracle/{capture.wl, capture.py, verify.jl}`

Extended for the new pole-bridge demo points `z ∈ {0.5, 0.95, 1.05,
1.4}` and the BF-256 string at `z = 1.05` (80-decimal-digit
`WeierstrassP[1.05 + c1, {0, c2}]` rendered via `N[..., 80]`).  Old
FW Table 5.1 v2 emissions are kept verbatim — the v1 pivot does not
invalidate them; they're now just out-of-scope for the v1 acceptance
and remain pinned for the v2 path-network work item
(`padetaylor-8cr`).

`mpmath.odefun` cleanly integrates `u'' = 6u^2` to `z = 0.95` at 40
dps; integration restarts from a fresh interpolant past the pole are
not needed for the v1 oracle (only `z < 1` covers the mpmath cross-
check).  Mathematica's `NDSolve` at WP=50 carries the cross-check
past the pole to `z = 1.4` with zero observable disagreement against
the closed-form (rel-err `0` to print precision — the `NDSolve`
restart at `z = 1.05 → z = 1.4` was clean, no convergence warnings).

### 2. `src/Problems.jl` — was 15-LOC stub; now 261 lines / ~116 effective LOC

The shipped public driver layer:

  - `struct PadeTaylorProblem{F, T, Y}` — wraps `f`, coerced `y0`,
    `zspan`, `order`.  `Y = Tuple{T, T}` selects the 2nd-order branch;
    `Y = T` is reserved for a 1st-order branch deferred to v2 (throws
    a `Suggestion` from `solve_pade` until then).
  - `struct PadeTaylorSolution{T, Y, P}` — trajectory + per-segment
    Padé store.  Callable as `(sol)(z) -> (u, u')` for the 2nd-order
    branch.
  - `solve_pade(prob; h_max, max_steps = 100_000)` — fixed-`h_max`
    driver consuming `pade_step_with_pade!` segment-by-segment;
    accumulates `(z[k], y[k], h[k], pade[k])`.
  - `taylor_eval(coefs, h)` — Horner-evaluator for the side-by-side
    pole-bridge demo (test 6.1.5).  ~5 LOC including bounds.

Fail-fast contract per CLAUDE.md Rule 1: throws with a `Suggestion`
on `order < 2`, degenerate `zspan`, `h_max ≤ 0`, evaluation outside
`[z_start, z_end]`, and exhausted `max_steps`.  Numerical breakdowns
inside the stepper (singular `C̃`, `Q(t) = 0` at evaluation point)
propagate from the lower layers unchanged.

The top-of-file docstring expands the v1/v2 acceptance scope
distinction and the `t = (z - z[k]) / h[k]; u' = P'(t) / h_k`
chain-rule formula for dense interpolation.  No content duplicated
from worklog 004 — the docstring points at the worklog and the
underlying `PadeStepper`-docstring rationale.

### 3. `src/PadeTaylor.jl` — added `taylor_eval` to umbrella re-exports

Single-line addition to the `using .Problems: …` line and the matching
`export` statement.  Without it the headline test 6.1.5 cannot
write `PadeTaylor.taylor_eval(coefs, 1.05)` cleanly.

### 4. `test/problems_test.jl` — replaced placeholder with six testsets

Testsets 6.1.1 – 6.1.6 implementing the test plan from
`docs/worklog/004-phase-6-pivot.md:112-122`.  Test 6.1.7 (mutation-
proof procedure) is documented inline in
`test/problems_test.jl:114-126` as a comment block matching the
Phase-5 pattern; the actual mutations were applied and verified by
the orchestrator before close-out (see "Mutation-proof procedure"
below).

## Phase-6 acceptance: 218 / 218 GREEN

Empirically observed relative errors against the three-source oracle:

| test | quantity | rel-err | rtol | margin |
|---|---|---|---|---|
| 6.1.1 (Padé, z=0.5)  | `u`  | 2.22e-16 | 1e-13 | ~3 orders |
| 6.1.1 (Padé, z=0.5)  | `u'` | 4.45e-16 | 1e-13 | ~3 orders |
| 6.1.3 (Padé, z=1.05) | `u`  | 3.45e-10 | 1e-9  | ~½ order |
| 6.1.3 (Padé, z=1.05) | `u'` | 4.52e-10 | 1e-9  | ~½ order |
| 6.1.5 (Padé,   z=1.05)        | `u` | 3.45e-10            | 1e-9             | ~½ order |
| 6.1.5 (plain Taylor, z=1.05)  | `u` | **2.50** (250%)     | rel-err > 0.1    | ~24× |
| 6.1.6 (BF-256, order=40, z=1.05) | `u` | 4.59e-15        | 1e-13            | ~2 orders |

**Headline.**  At `z = 1.05` (just past the pole at `z = 1`), Padé
beats plain Taylor by ≈ 9.86 orders of magnitude on identical Taylor
coefficients.  Tests 6.1.3 + 6.1.5 are the compelling pair: the same
input jet, one path through `robust_pade` and one through naïve
truncation, two answers separated by ten orders of magnitude.  This
is the analytic-continuation property of the algorithm in one test
setup.

## Algorithmic finding — the (15, 15) Padé approximation-error floor at t = 0.7

The test-spec drift caught at impl time, in the lineage of worklog
002's `0.7`-safety-factor and `/e²` corrections.

The original test-plan draft for 6.1.6 (per
`docs/worklog/004-phase-6-pivot.md:121`) set `order = 30` and
`rtol = 1e-13`, expecting BF-256 to clear Float64's `1e-9` result by
4 orders.  The subagent's GREEN run hit `3.46e-10` — a clean RED on
the BF-256 test, identical to Float64's number to plotting
precision.  Diagnosis (verified by an empirical orchestrator probe
sweep, not band-aided per CLAUDE.md Rule 2):

```
Float64 sweep at z = 1.05:
  order = 30 → 3.45e-10
  order = 35 → 3.45e-10  (Padé approx-error floor; not arithmetic noise)
  order = 40 → 1.63e-08  (DEGRADED — SVD ill-conditioning of the 41×41 block)
  order = 50 → 1.22e-04  (catastrophic)

BF-256 sweep at z = 1.05 (16-digit FW ICs):
  order = 30 → 3.46e-10  (same Padé floor as Float64 — BF-256 buys nothing here)
  order = 35 → 2.99e-12
  order = 40 → 4.59e-15  (saturated near IC-precision floor)
  order = 50 → 5.54e-15
```

Two readings of this table:

  - At `order = 30`, BF-256 adds **no value** over Float64.  Both
    saturate at the `(15, 15)` Padé approximation-error floor at
    `t = 0.7` (in the rescaled variable; the pole sits at
    `t ≈ 0.667`).  The floor is an algorithmic property of the rational
    approximant — bounded above by the truncation-error of the best
    `(15, 15)` rational approximation to `℘` on the disc of radius
    `0.7` minus `0.667` around `t = 0.7`.  It does not move with
    arithmetic precision.

  - At `order = 40`, BF-256 is **required**.  Float64 has degraded by
    ~3 orders due to SVD conditioning of the 41-by-41 Toeplitz block
    (a thousand-condition jump per added row of order); BF-256
    cleanly converges to ~5e-15, the IC-precision floor for the
    16-digit FW initial conditions after amplification through the
    pole.

The fix was to bump 6.1.6 to `order = 40` so the test demonstrates
BF-256's actual value-add over Float64 — the ratio `1.63e-8 / 4.59e-15
≈ 3.5e6` makes BF-256 a genuine accuracy enabler at this order, where
at `order = 30` it was a no-op.  Recorded as design rationale inline
in `test/problems_test.jl:91-100`.

The lesson is the worklog 002 lesson recurring: **specific numerical
spec values (rtol, order, tolerance bounds) drift from memory; only
empirical sweeps confirm what is actually achievable**.  Worklog 004
deliberately did not pin a rtol on 6.1.6, leaving the choice to the
implementer; that decision turned out to be the right call.  The
`(15, 15)` Padé approximation-error floor at `t = 0.7` is a real
algorithmic property of the rational approximant near the pole, not a
precision artefact, and it caps the achievable accuracy regardless
of the ground-arithmetic.

## Mutation-proof procedure (verified before close-out)

Both planned mutations from
`docs/worklog/004-phase-6-pivot.md:244-261` were applied; both
observed RED; both restored before the orchestrator's commit.

### Mutation A — drop the Padé conversion (use truncated Taylor)

`src/PadeStepper.jl:246`: replaced

    P_u = robust_pade(coefs_u_scaled, m, n)

with a stand-in that returns the trivial-denominator approximant
(numerator = `c̃`, denominator = `[one(T)]`, equivalent to truncated
Taylor on the rescaled variable):

    P_u = PadeApproximant{T}(coefs_u_scaled, [one(T)],
                             length(coefs_u_scaled) - 1, 0)

**RED bite**: `6.1.1, 6.1.2, 6.1.3, 6.1.4, 6.1.5 (Padé arm), 6.1.6` —
all six Phase-6 testsets bit at one or more `@test` calls.  Plus
**Phase 5** `5.1.1, 5.1.2, 5.1.3, 5.1.4` because `pade_step!` is the
underlying primitive there too.  The wide bite confirms `robust_pade`
is load-bearing on every Padé-Taylor step in the codebase, not just at
the demo `z = 1.5` segment.

### Mutation B — drop dense interpolation in `PadeTaylorSolution`

`src/Problems.jl` callable: replaced the 5-line dense-eval block

    h_k = sol.h[k]
    t   = (z_T - sol.z[k]) / h_k
    P_u = sol.pade[k]
    u   = _evaluate_pade(P_u, t)
    up  = _evaluate_pade_deriv(P_u, t) / h_k
    return (u, up)

with `return sol.y[k]` (segment-start state, ignoring the dense
interpolant).

**RED bite**: `6.1.1, 6.1.2, 6.1.3, 6.1.4, 6.1.5 (Padé arm only —
the Taylor arm reads `coefs_u` directly, bypassing `sol`), 6.1.6`.
**Phase-1 through Phase-5 untouched** (the mutation is local to
`Problems.jl`; lower layers do not depend on `PadeTaylorSolution`'s
callable).  Confirms the dense-interp formula `t = (z - z[k]) / h[k];
u' = P'(t) / h_k` is load-bearing — the test corpus genuinely
exercises dense interpolation, not just segment-endpoint evaluation.

## Frictions surfaced

### F1. The order/rtol coupling in 6.1.6 (described above)

Caught only by sweeping order across `{30, 35, 40, 50}` at both
Float64 and BF-256, not by reasoning from the spec.  The spec
deliberately under-pinned the test (worklog 004 did not fix a rtol);
the right next move was a sweep, not a spec amendment.  Leaving
implementer discretion on rtol-pinning was correct.

### F2. mpmath.odefun bottleneck was integration *range*, not mpmath itself

`mpmath.odefun` integrating `u'' = 6u^2` from `z = 0` to `z = 0.95`
at 40 dps takes ≈ 2 s — well within budget.  The previous
`capture.py`'s long hangs were caused by integrating to `z = 30` across
12 lattice poles, not by mpmath itself.  Confirmed by isolating the
short-range probe to `z ≤ 0.95`; the Phase-6 oracle infrastructure now
runs end-to-end in seconds.  No friction at the v1 scope.

### F3. UTF-8 mangling in Mathematica `WriteString` comment

The `∈` character (U+2208) in a source comment of `capture.wl` got
written as `a^\210\210` (octal escape of the UTF-8 bytes) into the
generated `oracles_wolfram.txt`.  Cosmetic only — affects only
comment lines in the generated text, no parser impact.  Lesson for
future probes: keep `WriteString` comment strings ASCII, or use
`Print` for diagnostic output and reserve `WriteString` for the
parsed payload only.

### F4. Cross-module private-function imports

`src/Problems.jl:88-89` imports `_evaluate_pade` and
`_evaluate_pade_deriv` from `PadeStepper`.  Underscore-prefixed names
are by convention "private but accessible"; this works in Julia
without ceremony.  Long-term, if a third caller appears (e.g. a
plotting helper, or a `recipe` for `Plots.jl`), these should be
promoted to public API in `PadeStepper`.  Recorded as an internal
note here rather than as a bead — premature abstraction risk per
CLAUDE.md Rule 9 ("v1-acceptable corner with the exact condition
that would force v2 work" = "third caller appears").

## Pointers

  - [`src/Problems.jl`](../../src/Problems.jl) — the shipped module +
    its literate top-of-file docstring.
  - [`src/PadeTaylor.jl`](../../src/PadeTaylor.jl) — umbrella; line 61
    + line 65 carry the new `taylor_eval` re-export.
  - [`test/problems_test.jl`](../../test/problems_test.jl) — six
    testsets; lines 91-100 carry the order/rtol-coupling rationale;
    lines 114-126 carry the inline 6.1.7 mutation-proof procedure.
  - [`external/probes/problems-oracle/`](../../external/probes/problems-oracle/)
    — extended capture (Mathematica closed-form + NDSolve, mpmath.
    odefun) + verifier; old FW Table 5.1 emissions kept verbatim for
    the v2 path-network work.
  - [`docs/worklog/001-stages-Z-1-2-handoff.md`](001-stages-Z-1-2-handoff.md)
    — original Stage-Z/1/2 handoff.
  - [`docs/worklog/002-phase-4-spec-correction.md`](002-phase-4-spec-correction.md)
    — the spec-drift pattern that 6.1.6's order/rtol coupling
    partially replicates.
  - [`docs/worklog/003-phase-5-padestepper.md`](003-phase-5-padestepper.md)
    — the underlying `pade_step_with_pade!` primitive that `solve_pade`
    composes; `h^k`-rescale and analytic-`u'` rationale.
  - [`docs/worklog/004-phase-6-pivot.md`](004-phase-6-pivot.md) —
    revised Phase-6 acceptance + test plan + oracle-extension
    checklist + subagent-prompt template that drove this shard.
  - `DESIGN.md §4 Phase 6` — revised acceptance.
  - `HANDOFF.md` — Phase-6 SHIPPED status block.

## Bead state

  - `padetaylor-tzm` — closed (Phase 6 v1 shipped).
  - `padetaylor-8cr` — remains open as the next P0 work item: FW 2011
    Table 5.1 long-range integration via the §3.1 path-network and/or
    a step-size selector that composes with Padé-bridge stepping.
  - Phase 7 (`CommonSolveAdapter` ext) and Phase 8
    (`PadeTaylorArblibExt`) become the next ready beads.  Phase 7 is
    optional per `DESIGN.md §4`; Phase 8 lifts the BF-256 test corpus
    onto the Arb element type once `Polynomials.roots` over `Arb` is
    settled (the open Phase-4 friction bead).
