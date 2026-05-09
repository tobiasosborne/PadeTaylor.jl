# HANDOFF.md — for the next agent picking up PadeTaylor.jl

> **Read this top to bottom before doing anything else.** Then read
> `CLAUDE.md` (the discipline rules). Then `RESEARCH.md` §1.1, §2,
> §3.3, §5, §7. Then `DESIGN.md` §4. Then this file again.
>
> If you skim and skip the references, you will repeat mistakes the
> previous session already paid for. The frictions surfaced are
> recorded in `docs/worklog/001-stages-Z-1-2-handoff.md`.

## TL;DR — where we are

**Stage 0 (research) and Stage 1 (design) are complete.** Stage 2
(implementation) is in progress — Phases Z, 1, 2, 3, 4, 5 of 9 are shipped:

| Phase | Module | Status | Tests |
|---|---|---|---|
| Z | scaffold + 3 ADRs + CLAUDE.md | ✅ shipped | 12 placeholders |
| 1 | `LinAlg.pade_svd` | ✅ shipped, mutation-proven | 14 |
| 2 | `RobustPade.robust_pade` | ✅ shipped, mutation-proven | 35 |
| 3 | `Coefficients.taylor_coefficients_*` | ✅ shipped, mutation-proven | 125 |
| 4 | `StepControl.step_jorba_zou` + `step_pade_root` | ✅ shipped, mutation-proven; DESIGN/HANDOFF spec correction | 7 |
| 5 | `PadeStepper.pade_step!` | ✅ shipped (the BIG inner-loop integration), mutation-proven | 16 |
| 6 | `Problems.PadeTaylorProblem`, `solve_pade` | ⏳ **NEXT — pivoted scope; see worklog 004** | spec'd in DESIGN.md §4 Phase 6 (revised) |
| 7 | `CommonSolveAdapter` ext | optional | DESIGN.md §4 Phase 7 |
| 8 | `PadeTaylorArblibExt` ext | optional | DESIGN.md §4 Phase 8 |
| 9 | Tier C: PI tritronquée pole-field qualitative | optional | DESIGN.md §4 Phase 9 |

**207 / 207 tests passing** as of the most recent Phase-6-prep commit.

**Phase 6 was pivoted** — the original FW 2011 Table 5.1 acceptance
(rel error `5e-13` at z=30 from FW ICs over a 60-segment integration
crossing 12 lattice poles) is **structurally unreachable** by the v1
algorithm.  FW's `7.62e-14` was achieved by their 5-direction
path-network (FW 2011 §3.1), which DESIGN.md §5 explicitly defers
to v2.  The revised Phase-6 v1 acceptance is a **compelling demo of
Padé-Taylor's analytic-continuation advantage over plain Taylor**,
verified against three-source consensus.  Full spec + test plan +
subagent-prompt template in `docs/worklog/004-phase-6-pivot.md`.

## Project shape

- `references/` — load-bearing PDFs (FW 2011, GGT 2013, FW 2014/15,
  RF 2014, FFW 2017, Jorba-Zou 2005, Mezzarobba 2019), each with
  `marker_single`-converted markdown under `references/markdown/<name>/`.
- `external/chebfun/padeapprox.m` — the canonical Padé reference impl.
  **Do not** modify; it is the oracle source.
- `external/probes/` — empirical probes that informed the design:
  - `taylorseries-arb/` — proves `TaylorSeries.jl::Taylor1{Arb}` works
    at order 80 (5/5 tests); `probe.jl` + `SYNTHESIS.md` are the
    record. **Phase 3 builds on this empirical confirmation.**
  - `padeapprox-oracle/` — `capture.m` runs `padeapprox.m` under
    Octave, writes `oracles.txt` (copied to `test/_oracles.jl`).
- `docs/adr/` — three accepted ADRs; read all three before code.
- `docs/worklog/` — frictions / lessons learnt; one shard so far.
- `src/` — six tier modules + the umbrella `PadeTaylor.jl`.
- `test/` — `runtests.jl` orchestrates per-module test files.
- `RESEARCH.md` — Stage 0 deliverable (1165 lines). Read §1.1 (FW 2011
  algorithm), §2 (GGT 2013 + reweighting), §3.3 (TaylorSeries+Arb
  empirical), §5 (bigfloat-SVD landscape), §7 (resolutions of PRD
  open questions). Skim §3 (adjacent packages) and §4 (verified
  arithmetic).
- `DESIGN.md` — Stage 1 deliverable. **§4 is the granular execution
  plan; treat as the source of truth for what each phase must build.**
- `CLAUDE.md` — 13 numbered discipline rules. Re-read after any
  context compression.

## Hard rules — do not deviate

1. **Ground truth before code.** Before writing a line of Phase 3 or
   later, open `RESEARCH.md` and the relevant FW 2011 / GGT 2013
   markdown line range. Cite by `references/markdown/<file>.md:<lines>`
   in commits.
2. **TDD always.** Tests RED before impl. For port-and-verify
   modules (Phase 4 step formulas), capture an oracle from the source
   (Octave for `padeapprox.m`; TaylorIntegration.jl for Jorba-Zou) and
   assert against it.
3. **Mutation-prove every load-bearing test.** Before committing a
   GREEN phase: perturb the impl in a *meaningful* way, confirm
   tests RED, restore. Document the procedure in test comments.
   The naïve mutation may not bite — see worklog shard 001 for the
   `GenericLinearAlgebra` piracy discovery.
4. **`bd create` any friction**. Beads is the only tracker.
   `bd ready` at start; `bd close <ids>` at end. Never `bd edit`.
5. **Commit after each GREEN phase** with a message that:
   - lists what was added (file paths + LOC)
   - lists tests added + which oracle was used
   - records the mutation procedure that bit
   - cites the references by line number
6. **Do not push to remote** — the user is doing local-first, no
   remote yet. Keep `git push` out of your Bash invocations.
7. **No parallel Julia agents.** One `julia` process at a time. The
   precompile cache is brittle to concurrent invocations.
8. **No author outreach.** Per project decision (CLAUDE.md Rule 12).

## Specific instructions for Phase 3 (`Coefficients`)

This is your immediate next task.

### Goal

Implement `PadeTaylor.Coefficients.taylor_coefficients_1st(f, z0, y0,
order)` and `taylor_coefficients_2nd(f, z0, y0, y1, order)` that
return the Taylor coefficients `c_0, …, c_order` of the local Taylor
expansion of the ODE solution at `z0`.

### Ground truth to read first

- `RESEARCH.md §1.1.1` — FW 2011 §2.1 method (b): "use one term,
  substitute, integrate, gain one correct coefficient per pass."
- `references/markdown/FW2011_painleve_methodology_JCP230/
  FW2011_painleve_methodology_JCP230.md:79-107` — the source.
  Specifically:
  - Lines 88-94 describe method (a) — symbolic differentiation. Don't
    use this for general `f`.
  - Lines 96-107 describe method (b) — operator-overloading. **Use this.**
- `RESEARCH.md §3.3` — empirical confirmation that
  `TaylorSeries.jl::Taylor1{T}` works for `T ∈ {Float64, BigFloat,
  Arb}`.
- `external/probes/taylorseries-arb/probe.jl` — concrete usage pattern.
- `external/TaylorIntegration.jl/src/integrator/jetcoeffs.jl` lines
  23-43 — the canonical coefficient-generation pattern in Julia for
  reference; we adapt the same pattern but expose it as a standalone
  function rather than as a step inside an integrator.
- `DESIGN.md §4 Phase 3` — the 5-test test plan.

### The algorithm (sketch — read FW 2011 §2.1.2 for full)

For first-order ODE `dy/dz = f(z, y)`:

1. Construct `z_taylor = z0 + Taylor1(T, order)` (a Taylor1 representing
   `z0 + h` symbolically, with `order` terms).
2. Bootstrap: start with `y_taylor = Taylor1([y0]; order)` (constant
   term only).
3. For `k = 1:order`:
   a. Compute `f_taylor = f(z_taylor, y_taylor)`. The result is a
      `Taylor1` whose `k`-th coefficient is now correctly computable
      from `y_taylor`'s first `k` coefficients (already known).
   b. Set `y_taylor[k] = f_taylor[k-1] / k` (the integration step:
      if `y' = f`, then `y_k = f_{k-1}/k`).
4. Return `y_taylor.coeffs` (the `Vector{T}` of length `order+1`).

For second-order ODE `d²u/dz² = f(z, u, u')`:

Similar bootstrap, but with two unknowns `u` and `up = u'` and the
relation `up_k = u_{k+1} · (k+1)` (formal differentiation of the
Taylor series). Update both per pass.

### Test plan (DESIGN.md §4 Phase 3)

| test | input | expected | shape |
|---|---|---|---|
| 3.1.1 | `f(z, y) = y`, `(z₀, y₀) = (0, 1)`, order 10 | `c[k] = 1/k!` (exp(z) Taylor) | spec-from-scratch |
| 3.1.2 | `f(z, y) = z² + y²`, `(z₀, y₀) = (0, 0)`, order 14 | matches `Taylor1` of the Bessel solution `t · J_{3/4}(t²/2) / J_{-1/4}(t²/2)` to ball-equality at order 14 | port-and-verify (FW 2011 §2.2.1) |
| 3.1.3 | `f(z, u, u') = 6u² + z` (PI), tritronquée IC, order 30 | finite, no NaN/Inf, leading coefficients agree with `WeierstrassP`-derived expansion to 1e-12 (sanity only — PI ≠ ℘ at large `z`) | spec-from-scratch |
| 3.1.4 | `T = BigFloat(precision=256)`, same as 3.1.1, order 60 | radii < `2⁻²⁰⁰` per `RESEARCH.md §3.3` empirical | port-and-verify |
| 3.1.5 | mutation-proof | replace one rule in TaylorSeries (e.g. negate `exp`) → 3.1.1 fails | RED |

For 3.1.2 the Bessel oracle:
- `J_{3/4}(t²/2)` at small `t` has Taylor expansion you can compute
  via `SpecialFunctions.besselj` over a `Taylor1`. Or use
  `Symbolics.jl` to derive the first 14 coefficients symbolically.
- Cleaner alternative: just integrate `dy/dz = z² + y²` symbolically
  by hand for the first ~6 coefficients (it's tractable: y_0=0,
  y_1=0, y_2=0, y_3=1/3, y_4=0, …).

### What `Coefficients.jl` should look like

Stub at `src/Coefficients.jl`:
```julia
module Coefficients
function taylor_coefficients_1st end
function taylor_coefficients_2nd end
end # module Coefficients
```

Replace with:
```julia
module Coefficients
using TaylorSeries: Taylor1
using ..PadeTaylor: ... (none needed at this layer)

export taylor_coefficients_1st, taylor_coefficients_2nd

function taylor_coefficients_1st(f, z0::T, y0::T, order::Int) where {T}
    # See module docstring for algorithm.
    z = z0 + Taylor1(T, order)
    y = Taylor1([y0], order)
    for k = 1:order
        f_t = f(z, y)
        y[k] = f_t[k-1] / k
    end
    return y.coeffs
end

function taylor_coefficients_2nd(f, z0::T, y0::T, y1::T, order::Int) where {T}
    # u'' = f(z, u, u'); two unknowns evolved jointly.
    # ...
end
end
```

**Watch the type stability.** `y0::T` should not be coerced to
`Float64` accidentally; use `T(y0)` explicitly when constructing
`Taylor1`. The empirical probe at `external/probes/taylorseries-arb/`
already validated this pattern works for `T = Arb`.

### What to commit

After Phase 3 GREEN + mutation-proven:
```
Phase 3 GREEN: Coefficients.taylor_coefficients_1st/_2nd; 5 tests + mutation-proof

src/Coefficients.jl (~80 LOC):
- ...

test/coefficients_test.jl:
- 3.1.1 ...
- 3.1.2 ...

Mutation-proven: ...

Test summary: 66/66 GREEN (61 + 5 new).
```

## Specific instructions for Phase 4 (`StepControl`)

### Status — Phase 4 SHIPPED 2026-05-09

Phase 4 is complete (commit recording StepControl, three-source oracle,
DESIGN.md correction). 192/192 tests GREEN. The original spec text
below contained an incorrect formula; the corrections (and the rest
of the paragraph) are kept here as a record of the spec drift caught
at implementation time. **Read `src/StepControl.jl`'s top docstring +
`docs/worklog/002-phase-4-spec-correction.md` before changing this
module.**

### Goal (as shipped)

Two step-size strategies in `src/StepControl.jl`:
- `step_jorba_zou(coefs, ε_abs; ε_rel = ε_abs)` — Jorba-Zou 2005
  **§3.3.1 eq. 11** in the fixed-order form (TI.jl-equivalent), NOT
  §3.2 eq. 3-8. The `0.7`-safety factor cited in the original spec
  has no source.
- `step_pade_root(P, z_current, target)` — FW 2011 §3.1 forward-
  projection-onto-direction heuristic.

### Ground truth (verified)

- `references/markdown/JorbaZou2005_taylor_IVP_package_ExpMath14/
  JorbaZou2005_taylor_IVP_package_ExpMath14.md:613-645` — §3.3.1 first
  step-size control, eq. 11.
- `external/TaylorIntegration.jl/src/integrator/stepsize.jl:17-89` —
  ported verbatim (`stepsize` + `_second_stepsize`).
- `external/probes/stepcontrol-oracle/{capture.jl, capture.py,
  capture.wl, verify.jl}` — three-source pin: TI.jl ≡ mpmath ≡
  wolframscript at 47 decimal digits on the canonical case.

### Polynomials.jl roots with Arb element type — still open

`Polynomials.roots` is validated for `Float64` / `Complex{Float64}`
only.  Arb element type for the Padé-root step is deferred to v2; see
the deferred bead `padetaylor-…` (file at session close per Rule 9).

## Phase 5 — SHIPPED 2026-05-09

Phase 5 (`PadeStepper.pade_step!`) is complete. Three-source oracle
(closed-form ℘ ≡ Mathematica `NDSolve` ≡ `mpmath.odefun`) at
<1e-15 relative agreement; mutation-proven on both `u` and `u'`
evaluation paths.  See `docs/worklog/003-phase-5-padestepper.md`
for the two algorithmic choices that crystallised at impl time
(`h^k` rescaling necessary near poles; analytic differentiation
of `P_u` for `u'` beats re-Padé by ~40×).

The Phase-5 `PadeStepper` now also exports
`pade_step_with_pade!(state, f, order, h) -> (state, P_u)` — used
by `solve_pade` for per-segment Padé storage.

## Specific instructions for Phase 6 (`Problems`) — REVISED 2026-05-09

**Read `docs/worklog/004-phase-6-pivot.md` top to bottom before
touching code.**  It has the full spec, test plan, oracle-extension
checklist, and a subagent-prompt template.

**The original FW 2011 Table 5.1 acceptance was abandoned.**  Two
subagents tried (fixed `h_max` + vault, and `step_jorba_zou`-adaptive
+ vault); both reach `0.45` relative error at z=30 vs the spec's
`5e-13`.  The remaining 12 orders of magnitude are algorithmically
inherent to fixed-real-axis stepping across multiple poles; FW's
published numbers used the 5-direction path-network deferred to v2.

**Revised v1 Phase-6 acceptance**: demonstrate that Padé-Taylor's
analytic continuation extends past the natural radius of convergence
of the underlying Taylor series — i.e. that one stored Padé bridges
a pole and gives correct values both before and *after* the pole.
Setup: build ONE Padé at z=0 with `h_max = 1.5`, spanning the ℘ pole
at z=1 (which sits at t=0.667 in the rescaled variable, strictly
inside the segment).  Evaluate `sol(z)` for z ∈ {0.5, 0.95, 1.05, 1.4}
and assert agreement with closed-form ℘ at progressively-loosened
tolerances.  Compare side-by-side with plain `taylor_eval` of the
underlying coefficients — Taylor diverges past z=1, Padé does not.

Three-source consensus:
- closed-form `WeierstrassP[z + c1, {0, c2}]` (Mathematica)
- `NDSolve` at WP=50 (Mathematica)
- `mpmath.odefun` at 40 dps (Python; only valid for z < 1)

Detailed test plan in worklog 004.  Oracle infrastructure already
partially in place (`external/probes/problems-oracle/`); needs
extension for the new test points (z = 0.5, 0.95, 1.05, 1.4).

**Deferred to v2** (file P0 bead at session close):
- FW 2011 Table 5.1 long-range integration. Requires the path-
  network or a step-size selector that adaptively shrinks h on
  the post-pole approach to keep `c̃`'s dynamic range tractable.

## Friction beads worth re-reading before starting

```
bd ready -n 30
```

The currently-open beads at session end:
- Phase 3-9 implementation tasks
- `padetaylor-<id>`: GenericLinearAlgebra pirates into LinearAlgebra.svd!
  (Stage 0 / Phase 1 finding)
- `padetaylor-<id>`: Polynomials.jl roots with Arb — Phase 4 probe
- Willers 1974 acquisition (low priority)

## Hard-won lessons (don't repeat these)

These come from the worklog shard `docs/worklog/001-stages-Z-1-2-handoff.md`.
Read that file in full before Phase 3.

Quick summary:
1. **`bd init` rejects dotted DB names.** `--prefix=padetaylor
   --database=padetaylor` is required for projects ending in `.jl`.
2. **`GenericLinearAlgebra.jl` pirates into `LinearAlgebra.svd!`.**
   Naïve `svd → LinearAlgebra.svd` mutation does not bite. Proper
   mutation for any "BigFloat dispatch" claim is `Matrix{Float64}(A)`
   downcast.
3. **Octave `i` ≠ Julia `im`** for the imaginary unit. Use
   `%.18eim` in `fprintf` format strings, not `%.18ei`.
4. **Mid-precision Padé table reduction is precision-dependent.**
   exp(z) (20,20) reduces to (7,7) at Float64 but stays at (20,20)
   at BigFloat-256. This is correct behaviour, not a bug.
5. **GGT 2013 §7 ill-posedness is real and visible empirically.**
   For log(1.2-z) (10,10), Octave-vs-Julia coefficients diverge at
   ~1e-3 *while the rational function values agree to 1e-15*. Test
   for function-value match, not per-coefficient match, in
   ill-posed-block cases.
6. **Test 1.1.1's mutation-proof procedure was wrong on first
   attempt.** Documented the *correct* mutation in the test comments.
   When you write Phase 3 mutation-proofs, verify the mutation
   actually bites before committing the procedure.

## Last commit before this handoff

`45e3e90` — Phase 2 GREEN. `git log --oneline` to see history.

Goodluck. Read CLAUDE.md again before you start.
