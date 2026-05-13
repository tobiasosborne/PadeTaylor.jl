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
(implementation) is in progress — Phases Z, 1–6 (original DESIGN.md
scope), plus the path-network and BVP modules (Phases 10–11, new
docs/design_section_6_path_network.md scope), are all shipped:

| Phase | Module | Status | Tests |
|---|---|---|---|
| Z | scaffold + 3 ADRs + CLAUDE.md | ✅ shipped | 12 placeholders |
| 1 | `LinAlg.pade_svd` | ✅ shipped, mutation-proven | 14 |
| 2 | `RobustPade.robust_pade` | ✅ shipped, mutation-proven | 35 |
| 3 | `Coefficients.taylor_coefficients_*` | ✅ shipped, mutation-proven | 125 |
| 4 | `StepControl.step_jorba_zou` + `step_pade_root` | ✅ shipped, mutation-proven | 7 |
| 5 | `PadeStepper.pade_step!` | ✅ shipped, mutation-proven | 16 |
| 6 | `Problems.PadeTaylorProblem`, `solve_pade` | ✅ shipped (Padé-vs-Taylor pole-bridge demo), mutation-proven; see worklog 005 | 11 |
| 7 | `CommonSolveAdapter` ext (`PadeTaylorAlg`) | ✅ shipped per worklog 009 | 34 |
| 8 | `PadeTaylorArblibExt` ext (SVD-only v1) | ✅ shipped per worklog 010 | 30 |
| **9** | Tier C: PI tritronquée pole-field qualitative | ✅ **shipped** — 4-of-5 pole-free sectors + conjugate symmetry + leading-pole magnitude at 25×25 [-4,4]² (worklog 012) | 12 |
| **10** | `PathNetwork.path_network_solve` | ✅ **shipped + FW Table 5.1 GREEN @ 2.13e-14 BF-256 (beats FW's 8.34e-14)**; see worklogs 006 + 008 | 33 |
| **11** | `BVP.bvp_solve` | ✅ **shipped (Tier-3 Chebyshev-Newton BVP)**, mutation-proven; see worklog 006 | 43 |
| **12 v1** | `Dispatcher.dispatch_solve` (1D chain) | ✅ **shipped (Tier-3 IVP↔BVP chain composition)**, mutation-proven; see worklog 007 | 45 |
| **12.5** | `EdgeDetector` — 5-point Laplacian pole-field classifier | ✅ **shipped** (worklog 011); FW 2011 §3.2.2 verbatim | 743 |
| **12 v2** | `LatticeDispatcher.lattice_dispatch_solve` — 2D-grid composition | ✅ **shipped** (worklog 013) — per-row BVP fill on smooth runs flanked by IVP cells (FW2011...md:190); FW Fig 4.1 quantitative pin deferred to follow-up `padetaylor-0c3` | 76 |
| 13 | `CoordTransforms` (FFW 2017 PIII/PV) | ⏸  P2, bead `padetaylor-bvh` | Tier-4 |
| 14 | `SheetTracker` (FFW 2017 PVI) | ⏸  P2, bead `padetaylor-grc` | Tier-5 |

**1377 / 1377 tests passing** as of the classical-Padé F64 default commit (worklog 021, bead `padetaylor-txg` closed).  **v0.1.0 tagged** (`38a49ae`, `CHANGELOG.md` ships); **Documenter site** generated (`30b3298`); **classical-Padé probe + ship** complete (worklogs 020 + 021) — F64 path now beats FW Table 5.1's z=10⁴ rel-err by 1.7× at full PathNetwork (1.4e-10 vs 2.34e-10); test acceptance tightened PN.2.2 z=30 F64 `1e-9 → 1e-12` and PN.2.3 z=10⁴ F64 `5e-5 → 5e-10`.

**Phase 6 shipped 2026-05-09 on the pivoted scope** — the v1
acceptance is a Padé-vs-Taylor pole-bridge demonstration (one stored
Padé bridges the lattice pole of `℘(z + c₁; 0, c₂)` at `z = 1`,
giving correct values at `z ∈ {0.5, 0.95, 1.05, 1.4}` while plain
truncation diverges past `z = 1`); see `docs/worklog/005-phase-6-
implementation.md`.  Pivot rationale + failure analysis in
`docs/worklog/004-phase-6-pivot.md`.

**Phases 10 + 11 shipped 2026-05-13** — Tier-2 path-network +
Tier-3 Chebyshev-Newton BVP solver, both GREEN + mutation-proven.
Three commits land (`910aab9` BVP ground truth; `0ada60f` Phase 10;
`cc7d8ca` Phase 11).  The session's deliverables include:

  - Phase 10 `PathNetwork.path_network_solve`: FW 2011 §3.1 5-direction
    wedge path-tree + Stage-2 fine-grid barycentric extrapolation.
    Generic in `T <: AbstractFloat` for `Float64` + `Complex{T}`.
    PN.1.1, PN.1.2, PN.2.1, PN.4.1 GREEN; PN.2.2 (FW Table 5.1
    quantitative ≤1e-13) + PN.3.1 (:steepest_descent test) deferred
    to follow-up bead `padetaylor-yt1`.
  - Phase 11 `BVP.bvp_solve`: Chebyshev spectral collocation per FW
    2011 §3.2 + Trefethen SMIM ch.6/13.  Generic in `T <: AbstractFloat`
    for `Float64` + `BigFloat-256` + `Complex{T}`.  Step-norm Newton
    convergence (eps^(3/4) default).  BV.1.1-BV.5.1 GREEN.
  - Full **figure-acceptance catalogue** at `docs/figure_catalogue.md`
    — 79 FW-family figures tiered T0-T5 with per-row acceptance.
  - **Unified path-network spec** at `docs/unified_path_network_spec.md`
    — 14 sections synthesised across all 5 FW-family papers.
  - **ADR-0004** at `docs/adr/0004-path-network-architecture.md`.
  - **BVP canonical recipe** at `references/bvp_recipe.md`
    — 8 sections covering Chebyshev nodes, D₁/D₂, Newton, barycentric.
  - **Octave oracle** at `external/probes/bvp-oracle/{capture.m, oracles.txt}`
    using DMSUITE chebdif/chebint; 47 pinned constants in
    `test/_oracle_bvp.jl`.
  - **Marker-converted primary refs**: Weideman-Reddy 2000 (DMSUITE),
    Trefethen SMIM 2000, Berrut-Trefethen 2004 (SIAM Review).  All
    paywalled, acquired via TIB VPN.
  - **DMSUITE source** at `external/DMSUITE/` (gitignored, MIT/public).

The session's worklog at `docs/worklog/006-phases-10-11-path-network-bvp.md`
covers orchestration pattern (Opus codes, Sonnet does grunt work),
algorithmic findings (step-norm vs residual-norm Newton convergence,
the Octave oracle's PI vs equianharmonic-℘ spec-drift catch), and
mutation-proof procedures for both modules.

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

## Phase 6 — SHIPPED 2026-05-09

Phase 6 (`Problems.PadeTaylorProblem`, `solve_pade`,
`PadeTaylorSolution`, `taylor_eval`) is complete on the pivoted v1
acceptance.  218/218 tests GREEN.  Three-source oracle (closed-form
`WeierstrassP` ≡ Mathematica `NDSolve` ≡ `mpmath.odefun`); BF-256
fourth source via `WeierstrassP[..., 80]` for the high-precision
test.  Mutation-proven on both the `robust_pade` call and the dense-
interpolation formula.

**Headline empirical result.**  At `z = 1.05` (just past the lattice
pole at `z = 1`), the Padé conversion of the local Taylor jet beats
plain Taylor truncation by ≈ 9.86 orders of magnitude on identical
input coefficients (rel-err `3.45e-10` vs `2.50`).  Tests 6.1.3 +
6.1.5 are the compelling pair — same jet, two paths, ten orders of
magnitude apart.

**Algorithmic finding caught at impl time.**  Test 6.1.6's BF-256
acceptance was bumped from `order = 30` to `order = 40`: at
`order = 30` both Float64 and BF-256 saturate at the same `(15, 15)`
Padé approximation-error floor (~3.5e-10 at `t = 0.7`), so BF-256
adds no value over Float64.  Only at `order = 40` does Float64
degrade to ~1.6e-8 (SVD ill-conditioning of the 41-by-41 Toeplitz
block) while BF-256 cleanly converges to ~5e-15 — making BF-256
required, not optional.  Full sweep + rationale in
`docs/worklog/005-phase-6-implementation.md` §"Algorithmic finding".

**Deferred to v2** (bead `padetaylor-8cr`):
- FW 2011 Table 5.1 long-range integration via the FW §3.1
  path-network and/or a step-size selector that composes with
  Padé-bridge stepping.  Failure analysis of fixed-`h_max` +
  pole-vault + Jorba-Zou-adaptive in
  `docs/worklog/004-phase-6-pivot.md` §"What was tried, what failed".

## Specific instructions for Phase 7 (`CommonSolveAdapter` ext)

This is your immediate next task.  Phase 7 is an optional package
extension wiring `PadeTaylorProblem` / `solve_pade` into the
`CommonSolve.jl` `init`/`step!`/`solve!` interface, so PadeTaylor.jl
plays cleanly with the wider `SciML` solver ecosystem.

**Read first** (in order):
1. `DESIGN.md §4 Phase 7` — Phase-7 acceptance + interface sketch.
2. `docs/adr/0003-*.md` — the ADR governing extension boundaries.
3. `src/Problems.jl` — the public driver this layer wraps.
4. `docs/worklog/005-phase-6-implementation.md` — the most recent
   shipped phase; orchestration pattern + mutation-proof discipline.

**Suggested workflow** (mirror Phase 3's "read first → test plan →
implementation skeleton" pattern from this HANDOFF lines ~96-225):

  - Stub `ext/CommonSolveAdapter/CommonSolveAdapter.jl` (a Julia
    weakdep extension on `CommonSolve`); update `Project.toml`'s
    `[weakdeps]` and `[extensions]` tables.
  - Sketch the test plan: at minimum, `init(prob; h_max) → integrator`,
    `step!(integrator) → (z, y)`, `solve!(integrator) → solution`,
    each cross-checked against the equivalent direct `solve_pade`
    call on the same problem.
  - Mutation-proof: drop the integrator's `pade` accumulator; assert
    the cross-check fails.

Phase 7 is **optional** per `DESIGN.md §4`: skip it if Phase 8
(`PadeTaylorArblibExt`) is the higher priority for downstream use.
The Phase-4 friction bead on `Polynomials.roots` over `Arb` element
type is a prerequisite for Phase 8 but unrelated to Phase 7.

## Friction beads worth re-reading before starting

```
bd ready -n 30
```

Open beads at end of PN.2.2-bugfix commit (11 total; `padetaylor-yt1`,
`padetaylor-1jf` closed):
- `padetaylor-8cr` **P0** — v2 umbrella FW Table 5.1 long-range
  (largely **subsumed** by PN.2.2 — u(z=30) now achieves 2.13e-14
  BF-256; bead may close after a confirmation pass at z=10⁴).
- `padetaylor-rgp` **P1** — Figure-acceptance catalogue (tracking).
- `padetaylor-c2p` **P1** — Edge detector (now consumed by Phase 12 v2 / `padetaylor-k31`).
- `padetaylor-k31` **P1** — Phase 12 v2: 2D lattice dispatcher with automatic edge detection.
- `padetaylor-kvi` **P1** — Phase 9 Tier C qualitative.
- `padetaylor-bvh` **P2** — Phase 13 CoordTransforms (Tier-4).
- `padetaylor-grc` **P2** — Phase 14 SheetTracker (Tier-5).
- `padetaylor-61j` **P2** — Willers 1974 acquisition.
- `padetaylor-8pi` **P2** — GenericLinearAlgebra svd! piracy friction.

**Next work item**: `padetaylor-8cr` (verify u(z=10⁴) = 21.0253033947…
matches FW Table 5.1 row 2 — should be straightforward now), then move
to Tier-3 figure reproduction via the Dispatcher.  Or attack figure
reproduction directly via `padetaylor-rgp` (e.g. FW Fig 4.6 BVP demo,
already in-tree as a BVP test).

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
7. **Empirically sweep `order` before pinning `rtol` on near-pole
   tests.** The `(m, n)`-Padé approximation-error floor depends on
   `m + n` and the distance to the next singularity, NOT on
   arithmetic precision; BF-256 only helps once you push past
   Float64's SVD-conditioning limit (typically `order ≥ 40`).
   Phase-6 test 6.1.6 hit this: at `order = 30` Float64 and BF-256
   both saturated at `~3.5e-10` (the `(15, 15)` Padé floor at
   `t = 0.7`); only at `order = 40` did BF-256 buy two extra
   orders of accuracy.  See `docs/worklog/005-phase-6-implementation
   .md` §"Algorithmic finding".

8. **Step-norm Newton convergence, NOT residual-norm, for spectral
   BVPs.**  The discrete residual `‖R‖_∞` floors at `cond(D₂) ·
   eps(T) ≈ N² · eps(T)` — an irreducible "spectral truncation"
   floor.  A residual-norm Newton tolerance below that floor is
   structurally unachievable; convergence loops forever / errors.
   The right choice is step-norm `‖Δu‖_∞ ≤ tol` with default
   `eps(T)^(3/4)` (production-Newton standard, NLsolve.jl-style).
   Phase-11 hit this on first cut; default was `100·eps(T) ≈ 2.2e-14`
   while PI N=20 floors at ~1e-12.  See
   `docs/worklog/006-phases-10-11-path-network-bvp.md` §"Step-norm
   Newton convergence".

9. **Oracle subagents can catch spec drift you didn't write.**
   Phase-11's Octave oracle subagent caught that the existing
   Phase-6 problems-oracle's `u_at_0_5 = 4.0044…` was for
   `u'' = 6u²` (equianharmonic ℘) NOT `u'' = 6u² + z` (Painlevé I).
   Recomputed via mpmath at 40 dps, cross-validated to ~1e-10.
   Lesson: oracle files should state their ODE prominently;
   `test/_oracle_problems.jl` line 1 should add a header comment
   distinguishing PI from equianharmonic ℘.  Recorded in worklog 006
   §"Spec-drift catch by the Octave oracle subagent".

10. **`marker_single` lives in `~/Projects/<project>/.venv/bin/`,
    not on system PATH.**  Same for many other Python tools (mpmath
    via uv, etc.).  Scan venvs before declaring a tool absent.  Saved
    as user memory `user-python-tools-venv-pattern`.

11. **TIB VPN gives access to paywalled SIAM/ACM/Elsevier journals.**
    Phase-11 acquired Weideman-Reddy 2000 (ACM TOMS) + Berrut-Trefethen
    2004 (SIAM Review) directly from publishers.  When user activates
    the VPN, attempt direct publisher URLs before preprint fallback.
    Saved as reference memory `reference-tib-vpn-paywalls`.

12. **ADR-specified invariants need invariant tests, not approximation
    tests**.  PN.2.2 surfaced a Phase-10 bug where the visited-node
    canonical-Padé invariant (specified in ADR-0004 + worklog 006) was
    NOT enforced in code — only the FW Table 5.1 long-range test at
    z=30 caught it (~22% rel-err Float64 AND BF-256).  No test asserted
    the invariant directly (e.g., "evaluating visited[k]'s stored Padé
    at t=0 returns visited_u[k]"); such a test would have caught the
    bug immediately at Phase-10 GREEN.  See worklog 008 §"What changed"
    and §"Hard-won lessons".

13. **When numerical error is invariant to precision, suspect the
    algorithm**.  PN.2.2's bug showed identical 22% rel-err at Float64
    AND BF-256.  Float64-rel-err == BF-256-rel-err is a strong signal
    that the issue is structural, not roundoff/truncation.  A
    pre-emptive BF-256 probe is worth the ~50× wall-time cost when a
    result smells wrong.  See worklog 008 §"Frictions surfaced" F3.

## Last commit before this handoff

Classical-Padé F64 default shipped (worklog 021, bead `padetaylor-txg` closed) — `src/RobustPade.jl` gains `classical_pade_diagonal` (~100 LOC) + element-type `method` dispatch.  F64 / Complex{F64} default `:classical`; BigFloat / Arb default `:svd`.  Singular Toeplitz auto-falls-back to `:svd`.  ADR-0005 documents.  Empirical: PN.2.3 z=10⁴ F64 rel-err `1.4e-10` (vs FW `2.34e-10` — beat by 1.7×) at full PathNetwork; worklog 020's `6.15e-11` was a wedge-walker-only number, Stage 2 interpolation adds ~2×.

### Session 2026-05-13 (post worklog 020) — `padetaylor-txg` closed (worklog 021)

  - **`src/RobustPade.jl`** — new `classical_pade_diagonal(c, m)` implementing FW 2011 §5.1.4 eqs. (5.4)+(5.5) via `lu(... ; check=false)` + `\`; throws `SingularException` on `!issuccess(F)`.  New `method::Symbol` kwarg on `robust_pade` with element-type-driven defaults (`_default_pade_method`).  Off-diagonal `(m≠n)` and singular Toeplitz both auto-route to existing SVD path.
  - **`test/classical_pade_test.jl`** (new, ~250 LOC) — 9 testsets / 63 assertions CP.1.1-9; mutation-proven (4 mutations P1-P4, all bite, all restored).  Bite cascades: P1 (sign-flip) cascades to 60+ RED across Phase 5/6/9 + PathNetwork + LatticeDispatcher; P3 (revert F64 default) bites 8 RED including the new rtol tightenings (proves they're load-bearing); P2 (drop fallback) bites 5; P4 (drop off-diagonal route) bites 4.
  - **`test/robustpade_test.jl`** — tests 2.1.2-4 + 2.1.6 (SVD-specific behaviour: diagonal-hop reduction, Froissart suppression, noise-thresholded recovery) explicitly pass `method=:svd` to keep testing GGT 2013 Algorithm 2 specifically.
  - **`test/pathnetwork_test.jl`** — PN.2.2 z=30 F64 rtol `1e-9 → 1e-12` (1000×); PN.2.3 z=10⁴ F64 rtol `5e-5 → 5e-10` (100,000×).  Both hold under classical.
  - **`docs/adr/0005-classical-pade-default-at-float64.md`** (new, ~150 lines) — element-type dispatch decision + when each method is right.
  - **`docs/worklog/021-classical-pade-default.md`** (new) — implementation log + mutation-proof bite counts + frictions.
  - Test suite: 1314 → **1377 GREEN** (+63 CP.1.* assertions).  Wall ~1m45s at 8 threads.

### Beads filed / closed this session

  - `padetaylor-txg` **closed** — classical-Padé F64 default shipped.

### Open beads end-of-session (worklog 021)

One P2 with concrete scope:
  - `padetaylor-7zw` P2 — BF-256 tritronquée pin (parallel testset for `test/fw_fig_41_test.jl`).

Plus existing P2/P3 admin / figure-pinning / perf beads (`padetaylor-bhw`, `-bho`, `-8ui`, `-1a3`, `-gky`, `-8pi`, `-61j`, and six P3s).  No P0 / P1 open.

### Hard-won lessons added this session (worklog 021)

**39. Mutation cascades expose architectural load-bearing-ness.** Mutation P1 (sign-flip on classical's RHS) was expected to bite ~10 CP.1.* assertions; it actually bit 60+ across Phase 5/6/9 + PathNetwork + LatticeDispatcher.  The mutation count is a direct measure of how deep an algorithm sits in the dependency graph; classical Padé is consumed by every step in every IVP path-network walk.  Lesson: don't pre-judge mutation bite scope — run it and see.

**40. Worklog 020's wedge-walker number is NOT the full-PathNetwork number.** Worklog 020 probed via a custom 5-direction wedge walker (Stage 1 only); the full `path_network_solve` adds Stage 2 barycentric interpolation, which is an extra Padé eval per target.  At z=10⁴ F64 this adds ~2× to the rel-err (6.15e-11 → 1.4e-10).  Lesson: when a probe simulates a stage, the production code's full pipeline carries additional error.  Always re-measure under the production code path before pinning test bounds.

**41. Degenerate test inputs are easier to construct by accident than expected.** CP.1.7's first attempt used `c_k = (1/2)^k` (geometric series) intending to demonstrate classical's no-reduction property — but the resulting Toeplitz is rank-1 (every row is a scaled copy of row 1), so classical correctly throws SingularException.  The "well-conditioned" claim of a Taylor series is a property of its convergence radius, NOT of the Toeplitz it builds.  Lesson: when picking a non-degenerate test case, compute or estimate the Toeplitz's rank explicitly before pinning expected behaviour.

### Session 2026-05-13 (post-v0.1.0, part 3) — classical-Padé Toeplitz `\` finding (worklog 020)

User-driven investigation: "did FW use arb prec?"  Answer: no.  FW
2011 line 350 explicitly describes their Padé method —
**Toeplitz approach + MATLAB `\` (LU/QR), no SVD anywhere**.  GGT
2013 (which we port via `src/RobustPade.jl`) was published two
years after FW 2011 and replaces `\` with SVD for robustness
against Froissart doublets.  For well-conditioned cases (℘ on the
equianharmonic trajectory), SVD's robustness is unused and costs
both per-step accuracy and wall time.

Empirical wedge-walker probe (5-direction `:min_u`, h=0.5, order=30,
F64):

| target | method     | wall (s) | rel-err   | vs FW Table 5.1   |
|--------|------------|----------|-----------|-------------------|
| z=30   | :svd       | 5.84     | 6.6e-12   | 87× worse         |
| z=30   | :classical | 0.01     | 1.54e-13  | 2× worse          |
| z=10⁴  | :svd       | 17.76    | 6.05e-6   | 25,800× worse     |
| z=10⁴  | :classical | 3.45     | 6.15e-11  | **3.8× BETTER**   |

Per-step probe (10-step trace along 22.5° wedge, F64 vs BF-256
truth): step 3 (closest to the z=1 lattice pole) — SVD 6.2e-10,
classical 7.1e-13 — **870× per-step accuracy improvement**.

**Beads filed (3)**:

  - `padetaylor-txg` **P1** — *Ship classical-Padé as F64 default*.
    Add `classical_pade_diagonal(c, m)` to `src/RobustPade.jl` per
    `FW2011_*.md:346-350`; dispatch by element type (classical for
    `Float64`/`Float32`/their complex, SVD-with-Jacobi for
    `BigFloat`/`Arb`); fallback to SVD on singular Toeplitz.
    Tightens PN.2.2 rtol `1e-9 → ~1e-12` and PN.2.3 rtol
    `5e-5 → ~1e-10`.  Acceptance + ADR-0005 sketched in the bead
    description.
  - `padetaylor-7zw` **P2** — *BF-256 tritronquée pin: rerun FW
    Fig 4.1 step-(i) BVP at BigFloat-256*.  BVP is already generic
    in `T <: AbstractFloat` (BV.5.1 confirms); test/fw_fig_41_test.jl
    just needs a parallel BF-256 testset.
  - `padetaylor-u7o` **superseded by padetaylor-txg** — the four
    speculative attack vectors are obsolete; classical-Padé is the
    answer.

No source changes committed yet — the implementation work is the
next session, blocking on the new P1 bead.

### Session 2026-05-13 (post-v0.1.0, part 2) — z=10⁴ investigation + PathNetwork enhancements (worklog 019, commit 29de074)

Pivot from naive "rerun PN.2.2 at z=10⁴ BF-256" (~4.5 h wall) to a
Float64 investigation question driven by what FW actually claimed.

  - **PN.2.3 testset** (`test/pathnetwork_test.jl`): routine Float64
    z=10⁴ regression test, `rtol = 5e-5` vs `u_at_10000_FW_ref =
    21.02530339471055`, plus `|imag(u)| < 5e-5` and `length(visited_z)
    > 20_000` (step-count regression detector).  Mutation E
    (canonical-Padé invariant) re-verified — bites PN.2.2 + PN.2.3 +
    cascade through Phase-9 / LatticeDispatcher.
  - **PathNetwork.path_network_solve enhancements**: new opt-in
    kwargs `verbose::Bool = false` + `progress_every::Integer = 500`
    for eager-flushed Stage-1 progress lines; `_wedge_evaluations`
    now runs the 5 wedge candidates under `Threads.@threads` (per-
    thread `PadeStepperState`, deterministic argmin by index).
    User RHS `f` must be thread-safe (documented in docstring).
    Measured 1.32× speedup at 8 threads on full `Pkg.test()`.
  - **F64 sweep over `(order, h, z)`** documented in worklog 019.
    Finding: `order` saturates by 30 in F64 long-range (truncation
    below roundoff); raising `order` does NOT close the FW gap.
  - **`external/probes/pathnetwork-long-range/capture.jl`**: offline
    BF-256 sweep probe over `z ∈ {30, 100, 500, 1000, 10000}`,
    ~3.3 h wall at 8 threads, re-runnable.  Not committed with full
    `oracles.txt` (multi-hour run; placeholder in commit).
  - `Project.toml`: add `Printf` (stdlib) for verbose-mode
    `@sprintf` / `@printf` helpers.

Test suite: 1311 → **1314 GREEN** (+3 PN.2.3 assertions).  Wall
~1m45s at 8 threads.

### Session 2026-05-13 (post-v0.1.0, part 1) — Documenter site + v0.1.0 release (commits 30b3298, 38a49ae, tag v0.1.0)

Two P1 beads closed:

  - **`padetaylor-36w`** — Documenter.jl docs site at `docs/{Project.toml,
    make.jl, src/*.md}`.  Self-bootstrapping build (`julia
    --project=docs docs/make.jl`); local-only per CLAUDE.md Rule 11
    (no `deploydocs`, no CI).  Sections: Home (overview + Phase-6
    pole-bridge headline + 14-tier status table), Architecture
    (synthesis of all four ADRs), API (per-module `@autodocs`
    blocks for all 14 modules + extensions), Figures (curated
    tier-by-tier from `docs/figure_catalogue.md`).
  - **`padetaylor-8ll`** — v0.1.0 release: Project.toml bumped
    `0.1.0-dev → 0.1.0`; new `CHANGELOG.md` with Keep-a-Changelog
    v0.1.0 entry covering all 14 phases + headline empirical
    results + known limitations; README.md status table refreshed
    (Phases Z-6 → all 14 tiers; tests 218 → 1311); annotated
    `git tag v0.1.0`.  **Not pushed** per CLAUDE.md Rule 6.

Local branch: 7 commits ahead of `origin/main`, tag `v0.1.0`,
nothing pushed.

### Beads filed this session (parts 1-3)

  - `padetaylor-36w`, `padetaylor-8ll`, `padetaylor-g9x` — all closed.
  - `padetaylor-txg` **P1** (classical-Padé F64 default) — open.
  - `padetaylor-7zw` **P2** (BF-256 tritronquée pin) — open.
  - `padetaylor-u7o` — superseded by `padetaylor-txg`.

### Open beads end-of-session (post worklog 020)

One P1 + several P2/P3:

  - **`padetaylor-txg`** P1 — classical-Padé default (next session).
  - `padetaylor-7zw` P2 — BF-256 tritronquée pin.
  - `padetaylor-g9x` (closed), `padetaylor-bhw`, `padetaylor-bho`,
    `padetaylor-8ui`, `padetaylor-1a3`, `padetaylor-gky`,
    `padetaylor-8pi`, `padetaylor-61j` — see `bd ready -n 30`.

### Hard-won lessons added this session (worklogs 019-020)

**32. Read the paper's wall times to constrain its arithmetic
precision.**  FW Table 5.1's 26.5 s for z=10⁴ Padé rules out
`vpa` — MATLAB's arbitrary-precision arithmetic is ~100× slower
than `double`.  Their numbers MUST be `Float64`.  When a paper
claims a specific rel-err number, the only-arithmetic-that-fits-
the-wall-time is often diagnostic of the underlying impl choices.

**33. Per-step roundoff amplifies non-linearly near pole crossings.**
10-step path-network probe at h=0.5 along 22.5° wedge sees two
algebraically-equivalent `h^k` rescaling schemes diverge by
2.5e-10 after just 10 steps — most concentrated in the single step
crossing the `z = 1` lattice pole (~65× amplification).  Intrinsic
to the flow near singularities, not a bug.

**34. `order` saturates at roundoff long before truncation in
Float64 long-range.**  Order 30, 40, 50 give identical rel-err at
z=10⁴ in Float64.  Raising `order` cannot close the FW gap; the
gap is roundoff-limited.

**35. When porting from a paper, port the actual algorithm, not
the "modern equivalent".**  We adopted GGT 2013 Algorithm 2 (Robust
Padé via SVD) — published two years AFTER FW 2011.  FW used
classical Padé via Toeplitz `\`.  Algorithmic progress moved the
recipe in a direction (SVD robustness) that costs accuracy on the
well-conditioned case.  When reproducing a paper's published
numbers, port the paper's algorithm.

**36. "Modern equivalent" sometimes means "added work the paper
didn't need".**  GGT 2013's robustness is real (Froissart
doublets, near-singular blocks).  On Painlevé pole-fields the
robustness machinery is unused; the bidiagonalization +
iterative diagonalization + tolerance-based degree-reduction cost
580–1000× per-step accuracy AND 5-580× wall time relative to
classical `\` on the smooth case.  Robustness has a price; spend
it only where needed.

**37. The right precision dispatch is element-type-driven.**  At
`Float64` the dynamic range is narrow and Demmel-Kahan SVD (or
LU) suffices.  At `BigFloat` it opens up and relative-accuracy
Jacobi SVD is the load-bearing tool.  Dispatch by `T`; default to
the right method per tier; expose alternatives as opt-in.

**38. Classical Padé via Toeplitz `\` is both faster AND more
accurate than GGT 2013 SVD for well-conditioned F64 Padé tables.**
Probe (z=10⁴, h=0.5, order=30, F64): SVD 17.8 s / 6.05e-6 rel-err;
classical 3.5 s / 6.15e-11 rel-err.  Classical beats FW's published
2.34e-10 by 3.8×.  See worklog 020 for full empirical table +
dispatch design (`padetaylor-txg`).

## Previous handoff content (pre-v0.1.0)

### Session 2026-05-13 late evening (part 4) — `padetaylor-grc` closed

One GREEN commit lands this session: Phase 14 / Tier-5 SheetTracker
for the sixth Painlevé equation per FFW 2017 §2.2.

  - `padetaylor-grc` **closed** — new module `src/SheetTracker.jl`
    (~150 LOC) exports:
    - `pVI_transformed_rhs(α, β, γ, δ) -> (ζ, w, wp) -> w''`
      Factory for the ζ-plane PVI RHS (FFW 2017 eq. 3, md:144).
      Removes the `z = 0` branch point; leaves `z = 1` as a
      `ζ = 2π·i·k` lattice on the imag axis.
    - `winding_delta`, `accumulate_winding`, `sheet_index` —
      path-side primitives for computing Riemann-sheet indices
      after a regular PathNetwork walk.  Normalisation to `(-π, π]`;
      sheet index via `round(total / 2π)`.
    - **PVI coordinate conversion = PV's** (FFW2017...md:146);
      callers reuse `pV_z_to_ζ` / `pV_ζ_to_z`.  No duplicate helper.
  - **26 new assertions** across 7 testsets (ST.1.1-ST.1.7):
    hand-pinned RHS at degenerate + non-degenerate parameters,
    end-to-end PVI direct-vs-transformed agreement at one Padé step
    (≤1e-10), winding primitives on closed CCW/CW loops + non-
    enclosing paths, sheet-index conversion, branch-lattice
    magnitude check at `ζ = 2π·im`.
  - Mutation-proven: Mutations O (RHS `/2 → /3`) + P (drop winding
    normalisation) + Q (`round` → `floor`) all bit as predicted.
  - `docs/figure_catalogue.md §6` row T5 marked PARTIAL: ζ-plane +
    winding primitives shipped; η-plane PVI eq. (Fig 2 column 1
    only) and constrained-wedge PathNetwork routing remain.

Test suite: 1285 → **1311 GREEN** (+25 ST.1.* + 1 umbrella `isdefined`).
Wall ~2m00s.

### Beads filed this session (part 4)

None.

### Open beads end-of-session (worklog 018)

Two P2 beads remain — both admin/non-code:
  - `padetaylor-61j` (Willers 1974 paper acquisition),
  - `padetaylor-8pi` (GLA piracy friction record).

No P0 / P1 / P2-code beads open.  All Tier-1 through Tier-5
architectural deliverables are at least PARTIAL acceptance.

### Hard-won lessons added this session (worklog 018)

29. **"Branch point in Float64" is a magnitude check, not isfinite**.
    At `ζ = 2π·im` exactly, `exp(ζ)` returns `1 + O(eps)·im` (the
    rounded `2π` constant is non-exact), so `(e^ζ - 1) = O(eps)` and
    the RHS blows up to `~1e30` — huge but still finite.  Test
    magnitude (`> 1e10`), not `isfinite`.  Downstream Padé-Taylor
    finite-Taylor-coefficient check is the actual fail-loud point.

30. **`atan2` angle differences need `(-π, π]` normalisation for
    winding**.  Raw `angle(z_new) - angle(z_old)` has `±2π` jumps
    at the branch cut.  Normalising per step to `(-π, π]` gives
    the "shortest signed path" interpretation, valid iff each path
    step is < π in angular extent.  Callers must walk with `h <<
    distance-to-branch` (worked out in the module docstring).

31. **Topological direction ≠ algorithm's same-side-wrap convention**.
    A step from `(-1, +0.01)` to `(-1, -0.01)` is geometrically
    downward (clockwise around origin) but normalised winding-delta
    reads `+0.02` (counterclockwise small).  The algorithm interprets
    every sub-π step as "going the short way around"; only actual
    circumambulation (loops > π built up over multiple steps)
    registers sheet changes.  Don't pin direction on individual
    cross-cut steps in tests.

## Last commit before this handoff (previous session)

CoordTransforms PIII/PV GREEN (bead `padetaylor-bvh` closed) — worklog 017.

### Session 2026-05-13 late evening (part 3) — `padetaylor-bvh` closed

One GREEN commit lands this session: Phase 13 / Tier-4 exponential
coordinate transforms for PIII and PV per FFW 2017 §2.1.

  - `padetaylor-bvh` **closed** — new module `src/CoordTransforms.jl`
    (~100 LOC) exports:
    - `pIII_transformed_rhs(α, β, γ, δ) -> (ζ, w, wp) -> w''` and
      same for PV — RHS factory closures ready for `PadeTaylorProblem`.
    - `pIII_z_to_ζ` / `pIII_ζ_to_z` and `pV_z_to_ζ` / `pV_ζ_to_z` —
      forward + inverse IC conversions at a single point.
    Module-header chapter cites FFW2017...md:39-48 verbatim (the two
    transformed equations P̃_III and P̃_V) and documents the
    deferrals (non-uniform Stage-1 nodes, adaptive Padé `h`, multi-
    sheet recovery).
  - **23 new assertions** across 6 testsets (CT.1.1-CT.1.6):
    IC round-trips, hand-pinned symbolic RHS values, end-to-end
    agreement direct-vs-transformed at one Padé step (≤1e-10).
  - Mutation-proven: Mutations L (PIII RHS `/4 → /3`) + M (PIII IC
    sign flip) + N (PV RHS `(w+1)/(w-1)` swap) all bit as predicted.
  - `docs/figure_catalogue.md §6` row T4 marked PARTIAL — ExpCoords
    shipped; non-uniform nodes + adaptive `h` are independent
    follow-ups, not filed as beads.

Test suite: 1262 → **1285 GREEN** (+22 CT.1.* + 1 umbrella `isdefined`).
Wall ~1m56s.

### Beads filed this session (part 3)

None.

### Open beads end-of-session (worklog 017)

Three P2 beads remain — same set as after worklog 016, minus `padetaylor-bvh`:
  - `padetaylor-grc` (Phase 14 SheetTracker — PVI, Tier-5),
  - `padetaylor-61j` (Willers 1974 acquisition),
  - `padetaylor-8pi` (GLA piracy friction).

No P0 or P1 beads open.

### Hard-won lessons added this session (worklog 017)

26. **Hand-derive IC-conversion factors on paper before coding**.
    PIII's `w' = (z u + z² u') / 2` has a load-bearing `z²` (not `z`)
    on the `u'` term, from the `dz/dζ = z/2` factor in the chain rule.
    Pin the result in a test with a hand-computed value (CT.1.1's
    `5/8`); mutation-test the sign + factor.

27. **End-to-end agreement is a stronger spec than a symbolic RHS pin**.
    Symbolic pins (CT.1.2, CT.1.5) catch transcription errors with
    crisp messages; direct-vs-transformed agreement (CT.1.3, CT.1.6)
    catches algebra errors that compensate at one sample point.  Ship
    both.

28. **Helpers-only module beats premature driver abstraction**.  Phase
    13's scope could have been a full `pIII_pV_solve(...)` driver.
    Choosing helpers-only (RHS factory + IC conversion) lets callers
    compose with `PadeTaylorProblem` + `path_network_solve` naturally,
    ~`100 LOC` source vs ~`250 LOC` for a driver wrapper, same
    downstream usability.

## Last commit before this handoff (previous session)

FW 2011 Fig 4.1 step-(i) BVP quantitative pin GREEN (bead `padetaylor-0c3` closed) — worklog 016.

### Session 2026-05-13 late evening (part 2) — `padetaylor-0c3` closed

One GREEN commit lands this session: reproduce the tritronquée FW
eq. 4.1 reference values via the FW Fig 4.1 step-(i) recipe.

  - `padetaylor-0c3` **closed** — `test/fw_fig_41_test.jl` (new ~110
    LOC) demonstrates the canonical FW Fig 4.1 step (i): Chebyshev-
    Newton BVP for `u'' = 6u² + z` on `[-20i, +20i]` with leading-term
    `u(z) = -√(-z/6)` Dirichlet BCs.  At N=240, `u(0)` pins to ≤3.5e-13
    and `u'(0)` to ≤5.3e-11 vs FW eq. 4.1 — both well under the bead's
    1e-10 spec.  6 testsets / 12 assertions; structural mutation-self-
    proof at FF.1.6 (the `+√` branch diverges).
  - **No `src/` changes**.  The recipe is the test itself; encapsulating
    it in a `fw_fig_41_axis_bvp(z_max; ...)` wrapper would be premature
    abstraction (no downstream caller, the bead is a pin not a feature).
  - `docs/figure_catalogue.md` row 4.1 marked PARTIAL with the step-(i)
    pin shipped + deferral notes for steps (ii) and (iii).

Test suite: 1250 → **1262 GREEN** (+12 FF.1.*).  Wall ~1m49s.

### Beads filed this session (part 2)

None.

### Open beads end-of-session (worklog 016)

Four P2 beads remain — same set as after worklog 015, minus
`padetaylor-0c3`:
  - `padetaylor-bvh` (Phase 13 CoordTransforms),
  - `padetaylor-grc` (Phase 14 SheetTracker),
  - `padetaylor-61j` (Willers 1974 acquisition),
  - `padetaylor-8pi` (GLA piracy friction).

No P0 or P1 beads open.

### Hard-won lessons added this session (worklog 016)

23. **The leading-term Dirichlet BC pins `u` exactly but not `u'`**.
    For BVP step-(i) recipes that derive `u'(endpoint)` from the BVP
    (no Neumann BC), the converged `u'` differs from the analytical
    leading-term derivative by an o(1) correction.  At `z = ±20i` for
    tritronquée: gap ≈ 2.5e-4.  Use a looser tolerance documenting
    the gap, or fetch an independent oracle for a tighter pin.

24. **N-tuning probe before pinning rtol on a fresh BVP**.  Different
    segments need different N.  `[-18, -14]` (Fig 4.6): N=20 suffices.
    `[-20i, +20i]` (Fig 4.1): N=240 for ≤1e-10.  Probe spectral
    convergence empirically (each +20 in N typically buys ~10× error
    reduction) BEFORE committing test's N parameter.

25. **`residual_inf` is not a convergence metric for spectral BVPs**.
    The BVP's step-norm Newton (worklog 006 lesson 8) converges in
    `‖Δu‖_∞ ≤ eps^(3/4)` while `‖R‖_∞` floats at `N²·eps(T)` or higher.
    Gating tests on `residual_inf ≤ TINY` is structurally wrong —
    Newton CAN converge with `‖R‖_∞ ~ 1e-8` at N=200 and still deliver
    `u(0)` to 1e-13.  Gate on per-eval error against the oracle, not
    on the residual.

## Last commit before this handoff (previous session)

PathNetwork `enforce_real_axis_symmetry` kwarg GREEN (bead `padetaylor-dtj` closed) — worklog 015.

### Session 2026-05-13 late evening — `padetaylor-dtj` closed

One GREEN commit lands this session: the library-level cure for the
y-asymmetry bug surfaced in worklog 014.

  - `padetaylor-dtj` **closed** — `PathNetwork.path_network_solve` gains
    opt-in `enforce_real_axis_symmetry::Bool` kwarg.  When `true`, walks
    only upper-half + on-axis canonical representatives (collapses ±0.0im
    signed zeros via `complex(real(z), abs(imag(z)))`) and mirrors lower-
    half outputs via `conj`.  Bit-exact `u(z̄) = ū(z)` invariant; roughly
    2× faster than full-plane walk.  Validates IC-on-real-axis precondition
    (throws `ArgumentError` if `Im(z₀) > tol`, `Im(u₀) > tol`, or
    `Im(u'₀) > tol`); trusts caller on real-coefficient `f`.  Default
    `false` preserves FW 2011 algorithm verbatim.
  - **4 new testsets / 13 assertions** PN.6.1-PN.6.4, mutation-proven
    (Mutations G + I bite as predicted; see worklog 015).
  - `examples/tritronquee_3d.jl` simplified — drops the user-space upper-
    half walk + manual mirror in favor of the kwarg.  Net -15 LOC user
    code, same PNG output.
  - `src/PathNetwork.jl` gains a new module-header chapter "Schwarz-
    reflection symmetry (opt-in)" documenting the algorithmic correctness
    preconditions and the design choices (signed-zero collapse, recursive
    call into the unkwarged path, on-axis output handling).

Test suite: 1237 → **1250 GREEN** (+12 PN.6.* + 1 extra in PN.6.2 for
`length(visited_z) > 1` walked-tree check).  Wall ~1m47s.

### Beads filed this session

None.

### Open beads end-of-session (worklog 015)

Five P2 beads remain — same set as before, minus `padetaylor-dtj`:
  - `padetaylor-bvh` (Phase 13 CoordTransforms),
  - `padetaylor-grc` (Phase 14 SheetTracker),
  - `padetaylor-61j` (Willers 1974 acquisition),
  - `padetaylor-8pi` (GLA piracy friction),
  - `padetaylor-0c3` (Fig 4.1 quantitative pin).

No P0 or P1 beads open.

### Hard-won lessons added this session (worklog 015)

20. **Julia `isequal(+0.0, -0.0) == false`** — and the same for
    `+0.0im` vs `-0.0im`.  Dict-keying on `Complex{Float64}` with
    on-axis cells is a footgun.  `complex(real(z), abs(imag(z)))`
    canonicalizes to `+0.0im` and eliminates the ambiguity.  Caught
    in PN.6.1's first GREEN attempt; documented in worklog 015 + the
    `_solve_with_schwarz_reflection` helper's docstring.

21. **Visited-tree-property tests need grids that force walking**.
    PN.6.2's first cut used a grid entirely within `h` of the IC ⇒
    Stage-1 took 0 steps ⇒ `visited_z == [IC]` regardless of mutation
    ⇒ Mutation I didn't bite the test.  Rewrote with off-axis targets
    > `h` from IC to force multi-step walks.  Lesson: when asserting
    properties of the visited tree, design the grid to grow it.

22. **Promote validated user-space workarounds to library API**.
    Worklog 014's `examples/tritronquee_3d.jl` carried the Schwarz-
    reflection workaround as 15 LOC of grid-partitioning user code.
    This session lifts it to an opt-in kwarg — the workaround was
    already exercised on the canonical PI tritronquée case, so the
    library version inherits the empirical validation.  Test plan
    just adds the bit-exact-symmetry invariants the user-space code
    couldn't easily assert.

### Session 2026-05-13 evening — all 4 P1 beads closed

Four GREEN commits land in this session, knocking out every P1 bead:

  - `9d1cc4c` **EdgeDetector** (Phase 12.5, `padetaylor-c2p` closed):
    FW 2011 §3.2.2 5-point Laplacian pole-field classifier.  60-LOC
    module, 743 assertions, 2 mutations.  Threshold on `log₁₀|Δu|`
    per FW Fig 3.3 (not bare `|Δu|`; bead description had a misread,
    documented in worklog 011).
  - `1c7056e` **Phase 9 PI tritronquée qualitative** (`padetaylor-kvi`
    closed): 5 testsets / 12 assertions verifying 4-of-5 pole-free
    sectors, conjugate symmetry, leading-pole magnitude at 25×25
    [-4,4]².  Mutation-proven (Mutation C bites 7/12 by flipping the
    ODE z-sign).  No source changes per DESIGN.md §4 "no new code".
  - `cc8c804` **Figure-catalogue refresh** (`padetaylor-rgp` closed as
    living document): `docs/figure_catalogue.md` §6 + §8 brought
    current; PathNetwork/EdgeDetector/BVP/Dispatcher v1 marked
    shipped; Fig 3.1 row marked PARTIAL with Phase 9 acceptance.
  - `4db29b1` **Phase 12 v2 LatticeDispatcher** (`padetaylor-k31`
    closed): 2D-grid composition machinery (~115 LOC) — PathNetwork +
    EdgeDetector partition + per-row BVP fill per FW2011...md:190.
    4 testsets / 76 assertions, mutation-proven (Mutations E + F).
    Tested against cosh closed form to ≤ 1e-10.

Test suite: 404 → **1237 GREEN** (+833 assertions across the four
shipped phases).  Wall ~1m45s.

### Beads filed this session

  - `padetaylor-0c3` **P2** Phase 12 v2.1: FW Fig 4.1 quantitative pin
    (`u(20i)` ≤ 1e-10) for `lattice_dispatch_solve`.  Different
    compositional pattern (vertical BVP + two outward pole fields)
    than v1's line-190 horizontal-row algorithm.  See worklog 013
    §"v1 scope decision".
  - `padetaylor-dtj` **P2** PathNetwork: add `enforce_real_axis_symmetry`
    kwarg.  `path_network_solve`'s `shuffle(rng, targets)` breaks
    conjugate symmetry for real-coeff/real-IC problems.  Workaround
    lives in `examples/tritronquee_3d.jl`; library-level kwarg is the
    proper fix.  See worklog 014.

### Open beads end-of-session

All remaining beads are P2:
  - `padetaylor-bvh` (Phase 13 CoordTransforms),
  - `padetaylor-grc` (Phase 14 SheetTracker),
  - `padetaylor-61j` (Willers 1974 acquisition),
  - `padetaylor-8pi` (GLA piracy friction),
  - `padetaylor-0c3` (Fig 4.1 quantitative pin, new this session).

No P0 or P1 beads remain open.

### Hard-won lessons added this session (worklogs 011-013)

14. **Trust the source paper over the bead description**.  `padetaylor-c2p`'s
    bead said `|∇²u|/h² > 0.001`; FW2011...md:208 verbatim is
    `log₁₀|Δu| > 0.001` (the `/h²` is inside `Δu` per eq. 3.3).  See
    worklog 011 §"Threshold interpretation".

15. **Probe orientation before writing test assertions on 2D-grid output**.
    For PI `u'' = 6u² + z`, the tritronquée's pole-bearing sector is
    centred on the *positive* real axis (angle 0), not negative.  The
    ASCII-dump-first heuristic caught this in seconds vs the
    theory-first approach that would have wasted 30+ min.  See worklog
    012 §"A correction to my own initial assumption".

16. **v2-bead acceptance criteria often blend "machinery" and "specific
    case"**.  `padetaylor-k31` fused the 2D-lattice composition layer
    with the FW Fig 4.1 acceptance pin.  These are different
    deliverables (different compositional patterns).  v1 ships the
    machinery; v2.1 ships the specific case.  See worklog 013 §"v1
    scope decision".

17. **`shuffle(rng, targets)` in PathNetwork breaks conjugate symmetry
    even when the underlying ODE preserves it**.  For real-coeff
    real-IC problems, the shuffle's random target order creates an
    asymmetric visited tree whose Stage-2 cascade can blow `|u|` apart
    by 4-5 orders of magnitude at conjugate-pair cells.  Fix: walk
    only upper-half + mirror via `u(z̄) = ū(z)`.  See worklog 014 §
    "Bug 1".

18. **Use even N for 2D Cartesian grids centred on a symmetry axis**.
    With odd N, the centre cell sits on the axis and gets walked
    along a special direction, producing path-dependent discontinuity
    from off-axis cells.  Even N avoids the special-case row.  See
    worklog 014 §"Bug 2".

19. **FW Fig 4.1 geometry ≠ line-190 horizontal-row algorithm**.
    FW Fig 4.1's vertical BVP with asymptotic BCs serves *open*
    smooth regions extending to infinity; line-190's horizontal-row
    BVP fill serves *interior* smooth runs bounded on both sides by
    pole fields.  Pure tritronquée has no bridgeable interior runs
    — the horizontal-row cure does not apply.  See worklog 014
    §"BVP-cure exploration".

### Examples directory shipped

`examples/` ships two working scripts + its own Project.toml/Manifest.toml:
  - `tritronquee_3d.jl` — 3D surface + 2D heatmap + EdgeDetector mask
    for PI tritronquée at 500×500 over `[-20, 20]²` in ~17 s.  Uses
    upper-half walk + conjugate mirror (Bug 1 workaround) + even N
    (Bug 2 workaround).  Reproduces the canonical FW Fig 3.1
    tritronquée pole-field qualitatively.
  - `tritronquee_bvp_compose.jl` — investigative reference for the
    "BVP-fill cures tritronquée artifact" claim (does not — 0 BVPs
    trigger on the pure tritronquée geometry).

PNGs are gitignored; regenerate via `julia --project=examples
examples/<name>.jl`.

Goodluck. Read CLAUDE.md again before you start.
