# Coverage snapshot — `src/` line coverage

Bead `padetaylor-krgy.7` (Tier-1 §1.4 of `03_hardening_methodology.md`).
Tool + runbook: `test/coverage/run_coverage.jl`, `test/coverage/README.md`.

**Coverage is a necessary-not-sufficient floor.** It finds lines no test
runs. It does **not** measure whether a covered line is *tested* — that
(test strength) is the mutation gate's job (`test/mutation/`). Read both.
**Do not chase 100%**: defensive guards, `verbose`/`show` blocks, and
arb-precision-only paths are legitimately cold at Float64.

## Inaugural snapshot (2026-06-09)

- Suite at capture: **9193 pass / 10 broken / 0 fail** (`Pkg.test(coverage=true)`).
- Repo HEAD: `0a8b3a0`.
- **OVERALL: 2543 / 2644 executable lines covered = 96.18%.**
- Coverage.jl run from a throwaway temp env (NOT a suite dep; absent from
  `Project.toml`); `src/*.cov` aggregated then deleted (gitignored, never
  committed).

Re-measure: `julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'` then
`julia test/coverage/run_coverage.jl` (aggregate-only). See the runbook.

## Per-file table (sorted ascending by coverage)

Only files below 100% are individually triaged below. The remaining **26**
files are at **100%** (`BranchTracker, Coefficients, CoordTransforms,
EdgeDetector, EdgeGatedSolve, Laplace2D, LatticeDispatcher, LinAlg,
NoumiYamada, NoumiYamadaSymmetry, PadeTaylor, PainleveHierarchy, PainleveNamed,
PoleField, SharedPade, SharedPadeCellB, SharedPadeDispatch, SheetTracker,
StepControl, VectorBVP, VectorCoefficients, VectorPathNetworkStage2,
VectorPoleField, VectorProblems, VectorStepper`). 100% ≠ tested — those files
are the mutation gate's domain.

| file | covered/total | % |
|---|---|---|
| Diagnostics.jl | 2/13 | 15.4 |
| VectorPathNetwork.jl | 86/108 | 79.6 |
| VectorStepControl.jl | 24/29 | 82.8 |
| PainleveSolution.jl | 50/56 | 89.3 |
| RobustPade.jl | 101/112 | 90.2 |
| PathNetwork.jl | 237/257 | 92.2 |
| Dispatcher.jl | 81/87 | 93.1 |
| PainleveClosedForm.jl | 71/76 | 93.4 |
| VectorWedgeStep.jl | 60/63 | 95.2 |
| SharedPadeDefect.jl | 52/54 | 96.3 |
| Problems.jl | 54/56 | 96.4 |
| IVPBVPHybrid.jl | 143/147 | 97.3 |
| Painleve.jl | 63/64 | 98.4 |
| Heun.jl | 83/84 | 98.8 |
| PadeStepper.jl | 93/94 | 98.9 |
| BVP.jl | 161/162 | 99.4 |

## Triage of the uncovered regions

Verdict legend — **(a)** genuine blind spot (→ propose a bead); **(b)**
intentionally-uncovered: defensive Rule-1 guard / `verbose`-print / `show`
display / unreachable / arb-prec-only.

### Diagnostics.jl — 15.4% (lines 209-219) — **(b)**
The entire uncovered block is `Base.show(::IO, ::MIME"text/plain",
::DiagnosticReport)` — the REPL pretty-printer. No test calls `display`/`show`
on a report. Pure display formatting; low value. The 15.4% headline is an
artifact of a tiny file (13 executable lines) being mostly one `show` method —
the numerical `quality_diagnose` machinery lives in the Diagnostics **ext**,
not here. *Optional* low-priority bead: a smoke `show` test (assert the
header/percent strings render) — but display-only, no numerical invariant.

### VectorPathNetwork.jl — 79.6% (432, 742-746, 847-849, 866-872, 892-895, 902-905, 979-1004) — **(b)**
Almost all `verbose`-mode `@printf`/`println` progress blocks and the two
verbose helpers `_verbose_target_start_v` / `_verbose_step_v` / `_fmt_v`
(979-1004) — never run because tests pass `verbose=false`. Line 432 is the
8-arg backward-compat constructor (pre-`visited_jets`). The `on_target_failure
=:skip` recovery path (catch → `push!(failed_targets)`, ~857-861) IS covered.
Display/back-compat only; **no blind spot.**

### VectorStepControl.jl — 82.8% (204-208, 214) — **MIXED: (a) + (b)**
The Jorba-Zou `_second_stepsize` **fallback scan** (204-208: `for j in
1:(p-2) … h2 = max(h2, (1/vnorm(c_j))^(1/j))`) — a real numerical fallback
that fires only when the primary estimate returns `Inf`; no Float64 test
reaches it. Line 209-213 is the zero-norm Rule-1 `throw` **(b)**. The scan at
204-208 is a **genuine blind spot (a)** — see proposed bead COV-1.

### PainleveSolution.jl — 89.3% (236, 239, 244, 246, 252, 265) — **(b)**
`_raw_summary`/`_raw_span` methods for `PathNetworkSolution` and the generic
fallback, plus `Base.show` branches (empty-params, IC-only span). Display
summary formatting; **no blind spot.**

### RobustPade.jl — 90.2% (154-155, 176-180, 272, 398, 468, 487) — **MIXED**
- 154-155 (`default_tol(Float64/Float32)`), 176-180 (`_default_pade_method`
  Float32/Complex variants) — trivial type-dispatch one-liners; the Float64
  numerical path is tested. Low value **(b)**.
- 468 — `error(...)` "all denominator coeffs below tol; should have been caught
  by r≡0 case" — Rule-1 defensive guard **(b)**.
- **272** — the trivial `(m,n)=(0,0)` case returns `r(z) ≡ c_0`; **(a)**.
- **398** — the `r ≡ 0` special case (`max|c| ≤ tol·‖c‖∞` → zero approximant,
  `padeapprox.m:69`); **(a)**.
- **487** — degenerate-`a` branch (`a = [0]` when all of `a` < `ts`,
  `padeapprox.m:134`); **(a)**.
  These three are real `padeapprox.m`-ported edge cases — see bead COV-2.

### PathNetwork.jl — 92.2% (484-487, 646-647, 653-657, 713-722, 951, 960-984) — **MIXED**
- 484-487 / 646-647 / 653-657 / 960-984 — `verbose` print blocks + verbose
  helpers **(b)**.
- 713-722 — `_attach_diagnostics` `diagnose=true`-without-extension `throw`
  (Rule-1 guard for the missing `PadeTaylorDiagnosticsExt`) **(b)**.
- **951** — `_evaluate_candidates` `catch → evals[k] = (z_cur, Inf, 0, nothing)`:
  the pole-landing recovery (a failed step at a candidate yields a
  worst-possible eval). A real numerical recovery branch; **(a)** — bead COV-3.

### Dispatcher.jl — 93.1% (124, 313-318, 368) — **(b) (mostly)**
124 = ComplexF64 convenience constructor; 313-318 = two `elseif` kwarg-combo
branches of the BVP-segment `bvp_solve` call (specific `initial_guess`/`tol`
combinations not exercised) — *kwarg combinatorics*, low value; 368 = unknown-
segment-type `throw` (Rule-1) **(b)**. Borderline: a test passing exactly one
of `{initial_guess, bvp_tol}` to a BVP segment would close 313-318 cheaply —
folded into bead COV-4 (low priority).

### PainleveClosedForm.jl — 93.4% (99, 150, 181, 229, 235) — **(b)**
Scattered single lines in closed-form `℘`/asymptotic branches and a `show`/
summary path; spot-checked as display + degenerate-parameter guards. No
numerical blind spot of note (closed-form values are cross-validated elsewhere).

### VectorWedgeStep.jl — 95.2% (216, 379, 441) — **(b)**
216 = `Base.showerror(::IO, ::VectorWalkError)` (display); 379 = `catch →
return -Inf` pole-landing recovery in the wedge disc-radius selector (mirrors
PathNetwork 951; the vector twin); 441 = analogous guard. Recovery/display.
The 379 recovery is the vector analogue of COV-3 — covered there conceptually;
**(b)** here (the scalar COV-3 test will exercise the shared pattern).

### SharedPadeDefect.jl — 96.3% (80-81) — **(b)**
The `while … abs(cc[end]) < 1e-30; pop!(cc)` trailing-near-zero strip loop
body — never executes because test `Q` vectors have no trailing near-zeros.
A defensive normalisation; benign. *Optional* low bead: feed a `Q` with a
trailing ~0 coefficient to exercise the strip.

### Problems.jl — 96.4% (130-131) — **(a) low**
The **scalar** (`!(y0 isa Tuple)`) `PadeTaylorProblem` constructor branch
(`y0_T = T(y0)`). Tests build problems with the 2nd-order Tuple `y0`; the
scalar/1st-order construction path is uncovered. Trivial to cover; bead COV-5
(low — the field assignment is simple, but Rule 5 wants a real assertion that
a 1st-order problem round-trips its `y0`).

### IVPBVPHybrid.jl — 97.3% (663, 673, 763-765) — **MIXED**
663 = `nothing` branch of the asymptotic-IC ternary; 673 = the no-
initial-guess `bvp_solve` dispatch (kwarg combinatorics) **(b)**. **763, 765**
= the `re ≤ first(slice_re)` / `re ≥ last(slice_re)` **boundary brackets** of
the BVP-slice evaluator — real boundary logic for evaluating outside the
interior slice range; **(a)** — bead COV-6 (low).

### Painleve.jl — 98.4% (154) — **(b)**
`_ident3(a,b,c) = (a,b,c)` — the identity coordinate map for no-transform
equations; a trivial one-liner reached only via a specific dispatch. Benign.

### Heun.jl — 98.8% (400-403) — **(a) low**
The `else` (non-real-params) branch of the Heun-equation solver dispatch — the
**complex-parameter** Heun path (with `branch_points`/`branch_cut_angles`).
Tests cover the real-params branch. A genuine blind spot for complex Heun
parameters; bead COV-7 (low-medium — needs a complex-param Heun fixture).

### PadeStepper.jl — 98.9% (593) — **(b)**
`_eval_poly_at_one` function header line flagged; the body is exercised via
`ffw_truncation_error`. Coverage line-attribution noise, not a real gap.

### BVP.jl — 99.4% (456) — **(b)**
Newton non-convergence `throw(ErrorException(...))` (max-iter exhausted) —
a Rule-1 fail-loud guard. Intentionally cold (the BVP fixtures converge).

## Summary

- **96.18% overall** — a healthy floor. 26 / 42 files at 100%.
- The bulk of uncovered lines are **intentionally cold (b)**: `verbose`-print
  blocks, `show`/display methods, and Rule-1 defensive `throw`/`error` guards.
  Chasing them to 100% would mean writing display-smoke tests of near-zero
  numerical value — explicitly discouraged (§1.4 "do not chase 100%").
- **7 genuine blind spots (a)** are proposed as beads COV-1…COV-7 below
  (mostly LOW priority; the numerically load-bearing ones are COV-2's
  `RobustPade` `padeapprox.m` edge cases and COV-1's Jorba-Zou vector fallback).
- Reminder: a 100% file is **not** a verified file. Cross-reference the
  mutation gate (`test/mutation/`) — e.g. `Coefficients.jl` and `RobustPade.jl`
  are 100% covered yet carry KNOWN/tracked mutation **survivors** under
  `padetaylor-98pe` (a covered-but-weak line). Coverage and mutation are
  orthogonal; neither alone is sufficient.

## Proposed beads (genuine blind spots — orchestrator triages; some may dup)

- **COV-1** `VectorStepControl.jl:204-208` — Jorba-Zou vector `_second_stepsize`
  fallback scan untested. *Test:* construct a vector jet whose primary step
  estimate is `Inf` (e.g. a jet with a vanishing leading derivative) and assert
  `vector_step_jorba_zou` returns the scan-derived `h2 = (1/‖c_j‖)^(1/j)` for
  the expected `j`. (LOW-MED; numerical fallback — Rule 5 needs a value
  assertion, not just non-throw.)
- **COV-2** `RobustPade.jl:272,398,487` — three `padeapprox.m`-ported degenerate
  branches: `(0,0)` trivial `r≡c_0` (272), `r≡0` all-coeffs-below-tol (398),
  degenerate-`a` `a=[0]` (487). *Test:* inputs that hit each — a constant
  series → (0,0); an all-below-tol series → r≡0; a series forcing
  trailing-`a`-near-zero → degenerate `a`. Assert the returned approximant
  matches `padeapprox.m` (cross-validate). (MED; numerical core.)
- **COV-3** `PathNetwork.jl:951` (+ vector twin `VectorWedgeStep.jl:379`) —
  the pole-landing `catch → eval = (Inf,…)` recovery in `_evaluate_candidates`.
  *Test:* drive a candidate node onto a known pole and assert it is ranked
  worst (loses every max comparison) rather than crashing the walk. (LOW-MED.)
- **COV-4** `Dispatcher.jl:313-318` — BVP-segment `bvp_solve` kwarg-combination
  branches (one of `{initial_guess, bvp_tol}` set). *Test:* a `dispatch_solve`
  with a BVPSegment passing exactly one of the two. (LOW; kwarg combinatorics.)
- **COV-5** `Problems.jl:130-131` — scalar (1st-order, non-Tuple `y0`)
  `PadeTaylorProblem` constructor. *Test:* build a 1st-order problem and assert
  it stores `y0_T == T(y0)` and round-trips through a solve. (LOW.)
- **COV-6** `IVPBVPHybrid.jl:763,765` — BVP-slice evaluator boundary brackets
  (`re ≤ first` / `re ≥ last` slice). *Test:* evaluate at a `re` below the
  first and above the last slice; assert it clamps to slice 1 / slice n. (LOW.)
- **COV-7** `Heun.jl:400-403` — complex-parameter Heun solver dispatch branch.
  *Test:* a Heun fixture with complex parameters; assert it routes through the
  `branch_points`/`branch_cut_angles` path and produces a correct value. (LOW-MED.)

> NOTE for the orchestrator: COV-2 and COV-3 touch the numerical core and may
> overlap existing `padeapprox.m`/walk-robustness beads (e.g. the `xhjw`/`jznu`/
> `fmf8` family or `padetaylor-98pe`); check before filing. The display/guard
> regions (Diagnostics `show`, BVP/Dispatcher/RobustPade `throw`s,
> `verbose` blocks) are deliberately **not** proposed as beads — they are
> intentionally-uncovered per §1.4.
