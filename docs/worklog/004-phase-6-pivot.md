# Worklog 004 — Phase 6 pivot: from "FW Table 5.1" to "Padé > Taylor demo"

**Date**: 2026-05-09
**Author**: Claude (orchestrating)
**Scope**: Phase 6 of Stage 2 — `Problems` module.  This shard
records the **abandoned** original Phase 6 acceptance criterion (FW
2011 Table 5.1 long-range integration to 5e-13) and the **adopted**
revised acceptance: a compelling demonstration of Padé-Taylor's
analytic-continuation advantage over plain Taylor, verified against
multi-source ground-truth consensus.

> **For the next agent picking up Phase 6**: read this worklog top
> to bottom, then `HANDOFF.md` Phase-6 section, then `DESIGN.md §4
> Phase 6` (revised).  The detailed test plan + algorithm sketch +
> pinned-oracle infrastructure are all in this shard.

## TL;DR

DESIGN.md's original Phase-6 acceptance was **structurally
unreachable** by the v1 algorithm: FW 2011 Table 5.1's `7.62e-14`
relative error at `z = 30` was achieved by FW's **5-direction
path-network** (FW 2011 §3.1 + §5.3 lines 352-368), which DESIGN.md
§5 explicitly defers to v2.  Two subagents tried fixed-`h_max` +
pole-vault + adaptive `step_jorba_zou` and reached `0.45` relative
error at z=30 — a finite trajectory that survives pole crossings
but accumulates ~12 orders of magnitude of error vs the spec.

The user's redirection: do not chase FW Table 5.1 in v1.  Instead,
**demonstrate the Padé-vs-Taylor distinction** — that the Padé
representation extends past the natural radius of convergence of
the underlying Taylor series — on a single representative ODE,
verified against three-source closed-form / NDSolve / mpmath.odefun
consensus.  That is the v1 deliverable.  Long-range integration
(FW Table 5.1 reproduction) becomes a v2 P0 with the path-network.

## What was tried, what failed, why

### Attempt A — fixed `h_max = 0.5` (the spec)

With FW ICs from z=0 and `h_max = 0.5`, step 2 lands at z=1.0,
which is **exactly** the lattice origin pole of `℘(z + c₁; 0, c₂)`
on the (`c₁ = -1`, `c₂ = 2`) lattice.  The Padé denominator
correctly captures the pole as a denominator zero at `t = 1` in the
rescaled variable; evaluating at `t = 1` yields ~7e12 garbage.
Subsequent Taylor jets at the garbage IC have wildly skewed
coefficients; `RobustPade._trim_and_normalise`'s threshold
(`tol·‖c̃‖₂`, faithfully ported from Chebfun's `padeapprox.m:134`)
is dominated by the high-order coefficients and trims the whole
numerator to zero.  Trajectory cascades to `u ≡ 0` from step 3
onward.

### Attempt B — vault when denominator vanishes at `t = 1`

Detect `|Q(1)| < tol · ‖b‖₁`, multiply `h` by `vault_factor = 1.2`,
rebuild Padé.  Step 2 vaults to `h = 0.6` and lands at z = 1.1
(0.1 past the pole at z = 1).  `u(z=1.1) ≈ 100` correctly per
℘'s `1/(z-1)²` near-pole behaviour.  But step 3 from z = 1.1 with
`h_max = 0.5` has radius of convergence `≈ 0.1` (distance to the
pole at z = 1, which is now behind us).  `c̃ = h^k · c_k` spans 22
orders of magnitude; `RobustPade`'s trim collapses the numerator.
Cascade to `u ≡ 0` at step 3.

### Attempt C — `step_jorba_zou` for adaptive shrinkage

`step_jorba_zou(coefs, eps)` returns the bound `(eps/|c[k]|)^(1/k)`
which shrinks as `|c[k]|` grows near a pole.  Per `_take_one_step!`
the trajectory takes h ∈ {0.34, 0.22, 0.15, 0.10, 0.06, …},
asymptotically approaching z = 1 but never crossing.  After
`max_steps`, the loop exits without reaching the target.

The Jorba-Zou bound is appropriate for **smooth-ODE accuracy
control**, not for FW's Padé-bridge use case where the bigger step
is exactly what bridges the pole.  Mixing the two paradigms
fights the FW algorithm.

### Attempt D — `tol = 0` in `robust_pade` to disable trim

Disabling the trim threshold lets the trajectory survive past
multiple pole crossings.  `u(30) ≈ 0.597` vs expected
`1.095098255959744`; relative error `0.45`.  Finite, qualitatively
℘-like (matches the elliptic structure with peaks near each
pole), but accumulating ~12 orders of magnitude vs the spec's
`5e-13` budget.

The accumulating error is **algorithmically inherent** to fixed-h
real-axis stepping across multiple poles.  Each pole-vault loses a
few digits to the Padé extrapolation past `t = 1`; FW's
path-network avoids this by complex-plane detours.

## The revised Phase-6 acceptance

**Goal**: demonstrate that PadeTaylor.jl's Padé-Taylor pipeline
extends the domain of validity of the local Taylor jet past the
natural radius of convergence — i.e. that Padé bridges poles where
plain Taylor diverges.  Verified against three-source consensus.

**Setup** (single-segment, no compounding error):

```julia
fW(z, u, up) = 6 * u^2

# ONE step, h_max = 1.5, spans the lattice pole at z = 1.
prob = PadeTaylorProblem(fW, (u_0_FW, up_0_FW), (0.0, 1.5); order = 30)
sol  = solve_pade(prob; h_max = 1.5)

# sol.pade[1] is the Padé built at z = 0 with rescaling h = 1.5.
# The pole at z = 1 sits at t = 1/1.5 = 0.667 — strictly INSIDE
# the segment, NOT at the endpoint.  No vault needed.  sol(z) for
# z ∈ [0, 1.5] evaluates the SAME Padé at t = z / 1.5 ∈ [0, 1].
```

**Test plan** (rewrite of `test/problems_test.jl`):

| test | what it shows | oracle source |
|---|---|---|
| 6.1.1 | `sol(0.5)` ≡ ℘(0.5) to 1e-13 | closed-form ≡ NDSolve ≡ mpmath.odefun |
| 6.1.2 | `sol(0.95)` ≡ ℘(0.95) to 1e-9 (near pole) | closed-form ≡ NDSolve (mpmath stops near pole) |
| 6.1.3 | `sol(1.05)` ≡ ℘(1.05) to 1e-9 — **Padé continues PAST the pole** | closed-form ≡ NDSolve (mpmath fails past pole) |
| 6.1.4 | `sol(1.4)` ≡ ℘(1.4) to 1e-7 — Padé still tracks ℘ further out | closed-form ≡ NDSolve |
| 6.1.5 | Plain `taylor_eval(coefs_at_0, h=1.05)` is GARBAGE (rel err > 0.1) — Padé wins | side-by-side proof |
| 6.1.6 | BigFloat-256 version of 6.1.1 / 6.1.3 | closed-form (`WeierstrassP` at WP=80) |
| 6.1.7 | Mutation-proof — see procedure at end | RED |

**The compelling line**: tests 6.1.3 + 6.1.5 are the headline.  Same
Taylor coefficients; one path (Padé) gives the right answer past
the pole; the other path (plain truncation) diverges.

**Acceptance**: all 6 tests pass.  Total suite target ≈
207 + ~14 individual @test calls = 221+/221+.

**Deferred to v2** (file a P0 bead at session close):
- FW 2011 Table 5.1 long-range integration (z ∈ {30, 10⁴, 28.261}).
  Requires the path-network (FW 2011 §3.1) or, alternatively, a
  step-size selector that adaptively shrinks `h` on the
  *post-pole* approach to keep `c̃`'s dynamic range under control.
- `solve_pade` with `step_strategy = :jorba_zou` or `:pade_root` —
  the Phase-4 step controllers exist but DO NOT compose cleanly
  with FW's Padé-bridge.  v2 work item.

## What's already in place for the next agent

### `src/PadeStepper.jl` — `pade_step_with_pade!` exported

The first Phase-6 subagent added:

```julia
function pade_step_with_pade!(state, f, order::Int, h::Real) -> (state, P_u)
```

This is the underlying primitive `solve_pade` should call.  It does
the Phase-5 5-step algorithm (Taylor jet → h^k rescale → Padé →
evaluate at t=1 → mutate state) AND returns the Padé.
`pade_step!` is now a thin wrapper that discards the Padé.  Phase-5
tests stay GREEN.

### `external/probes/problems-oracle/`

Already pinned for the FW Table 5.1 z values (z ∈ {30, 10⁴,
28.261, 7.123}).  **Extend** for the new test points:

  - z = 0.5, 0.95, 1.05, 1.4 (the new pole-bridge demo).
  - The BigFloat-256 case: `WeierstrassP[1.05 + c1, {0, c2}]` at
    `N[..., 80]`.

`capture.wl` only needs 4 new lines; rerun and verify.

### `test/_oracle_problems.jl` — auto-generated

Re-running `external/probes/problems-oracle/verify.jl` after
extending `capture.wl` will refresh the pinned file.  Add the new
test points to the `verify.jl` cross-check assertions too.

### Three-source consensus discipline

Phase 6's two-source (closed-form ≡ paper FW reference) is
sufficient for the deferred FW Table 5.1 oracle.  For the new
pole-bridge tests, three-source consensus is achievable for
`z < 1` (closed-form ≡ NDSolve ≡ mpmath.odefun) and two-source for
`z > 1` (mpmath.odefun cannot integrate through poles reliably).
Document the asymmetry in `verify.jl`.

## Detailed orchestration plan for the next agent

**Step 1 — extend the oracle.**
  - Edit `external/probes/problems-oracle/capture.wl`:
    - Add `WeierstrassP` evaluations at z ∈ {0.5, 0.95, 1.05, 1.4}.
    - Add an `N[..., 80]` string at z=1.05 for BigFloat-256.
  - Edit `external/probes/problems-oracle/capture.py`:
    - `mpmath.odefun` integrate `u'' = 6u^2` from z=0 to z=0.95
      (stops short of the pole at z=1; do NOT integrate past).
    - Output values at z = 0.5, 0.95.
  - Run `wolframscript -file capture.wl` and `python3 capture.py`.
  - Edit `external/probes/problems-oracle/verify.jl` to assert
    cross-source agreement on the new test points.
  - Run `julia --project=. external/probes/problems-oracle/verify.jl`.
    Should print all-agreement messages and emit the updated
    `test/_oracle_problems.jl`.

**Step 2 — write the new tests.**
  - Replace `test/problems_test.jl` with the test plan from this
    worklog (6.1.1 - 6.1.7).
  - Verify tests are RED before implementing (`MethodError: no
    method matching PadeTaylorProblem`).

**Step 3 — dispatch a subagent to implement `src/Problems.jl`.**
  - Subagent prompt template lives below; copy + adapt.

**Step 4 — verify GREEN, mutation-proof.**

**Step 5 — commit + update HANDOFF.md + close `padetaylor-tzm`.**

**Step 6 — file a P0 bead for v2 FW-Table-5.1 work.**

## Subagent prompt template for Step 3

> You are implementing the revised Phase 6 of PadeTaylor.jl.  The
> repo is at `/home/tobiasosborne/Projects/PadeTaylor.jl`.  Read
> `CLAUDE.md`, `docs/worklog/004-phase-6-pivot.md` (this file), and
> `test/problems_test.jl` before touching code.
>
> Build `src/Problems.jl` with:
>   - `struct PadeTaylorProblem{F, T, Y}` carrying `f`, `y0`, `zspan`,
>     `order` (default 30).  `Y = Tuple{T, T}` for 2nd-order, `Y = T`
>     for 1st-order; autodetect from `typeof(y0)`.
>   - `solve_pade(prob; h_max, max_steps=100_000) -> PadeTaylorSolution`
>     using `PadeStepper.pade_step_with_pade!` per step.  For Phase 6
>     v1 demo: pure fixed-`h_max` stepping, no vault, no Jorba-Zou.
>     The test problem's `h_max = 1.5` puts the pole strictly INSIDE
>     the segment (case (b) in PadeStepper's docstring), so vault is
>     unnecessary.
>   - `struct PadeTaylorSolution{T, Y}` with `z, y, h, pade` vectors.
>     Callable: `sol(z) -> (u, up)` for 2nd-order (or `T` for 1st).
>     Dense interp: find segment `k`, evaluate `sol.pade[k]` at
>     `t = (z - z[k]) / h[k]`; `u' = P'(t) / h_k`.
>   - Export a small `taylor_eval(coefs, h)` helper for the
>     side-by-side comparison test 6.1.5.  ~5 LOC.
>
> Tests must pass: 207 prior + 6 testsets ≈ 221+/221+ GREEN.
>
> ≤ 200 LOC.  Literate top-of-file docstring documenting the
> Padé-vs-Taylor demonstration goal + the v1/v2 acceptance scope
> distinction (cite worklog 004).

## Mutation-proof procedure for the next agent

Each load-bearing path needs an independent mutation that REDs at
least one test:

  - Mutation A — drop the Padé conversion, use truncated Taylor:
    in `solve_pade` (or `pade_step_with_pade!`), replace the
    `robust_pade(c̃, m, n)` call with one that returns `r = c̃` (a
    Padé equivalent to truncated Taylor).
    Expected: tests 6.1.3, 6.1.4 fail (Taylor diverges past pole).

  - Mutation B — drop dense interpolation:
    in `(sol::PadeTaylorSolution)(z)`, return only the
    segment-start state `sol.y[k]` ignoring the Padé eval.
    Expected: 6.1.1 fails (z=0.5 returns u(0)=1.07 not u(0.5)=4.0).

Each mutation must be VERIFIED to bite (per worklog 001 friction)
before being recorded in the commit message.

## Pointers

- `docs/worklog/001-stages-Z-1-2-handoff.md` — original handoff.
- `docs/worklog/002-phase-4-spec-correction.md` — StepControl spec
  drift caught at impl time (similar pattern: paper line range
  triangulation ).
- `docs/worklog/003-phase-5-padestepper.md` — PadeStepper inner-
  loop integration; established the `h^k` rescaling rationale and
  the analytic-differentiation choice for `u'`.
- `src/PadeStepper.jl` — Phase 5 inner loop.  `pade_step_with_pade!`
  is the primitive Phase 6 should call.
- `external/probes/problems-oracle/{capture.wl, capture.py,
  verify.jl}` — oracle infrastructure, partially populated.
- `test/_oracle_problems.jl` — pinned values; needs extension for
  the new test points.
- `DESIGN.md §4 Phase 6` — revised acceptance.
- `HANDOFF.md` — updated Phase-6 status.

## Bead state

`padetaylor-tzm` reopened with revised description matching this
worklog's pivot.  v2 P0 bead filed for FW Table 5.1 long-range work.
